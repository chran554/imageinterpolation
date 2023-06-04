package main

import (
	"fmt"
	img "image"
	"image/color"
	_ "image/jpeg"
	_ "image/png"
	"math"
)

type ImageInterpolator struct {
	kernelX0          int
	kernelY0          int
	kernel            [][]color.Color
	interpolationFunc func(float64, float64, [][]color.Color) color.Color
}

func NewBicubicImageInterpolator() ImageInterpolator {
	return ImageInterpolator{
		kernelX0:          1,
		kernelY0:          1,
		kernel:            [][]color.Color{{nil, nil, nil, nil}, {nil, nil, nil, nil}, {nil, nil, nil, nil}, {nil, nil, nil, nil}},
		interpolationFunc: bicubicInterpolateColor,
	}
}

func NewBilinearImageInterpolator() ImageInterpolator {
	return ImageInterpolator{
		kernelX0:          0,
		kernelY0:          0,
		kernel:            [][]color.Color{{nil, nil}, {nil, nil}},
		interpolationFunc: bilinearInterpolateColor,
	}
}

func NewNearestNeighbourInterpolator() ImageInterpolator {
	return ImageInterpolator{
		kernelX0:          0,
		kernelY0:          0,
		kernel:            [][]color.Color{{nil, nil}, {nil, nil}},
		interpolationFunc: nearestNeighbourInterpolateColor,
	}
}

func main() {
	fmt.Println("Up and running...")
	upscale := 50

	// original, err := readImage("resources/original.png")
	original, err := readImage("resources/wiki.png")
	if err != nil {
		panic(err)
	}

	width := original.Bounds().Dx() * upscale
	height := original.Bounds().Dy() * upscale

	bilinearInterpolatedImage := img.NewRGBA(img.Rect(0, 0, width, height))
	nearestNeighbourInterpolatedImage := img.NewRGBA(img.Rect(0, 0, width, height))
	bicubicInterpolatedImage := img.NewRGBA(img.Rect(0, 0, width, height))

	InterpolateImage(original, bicubicInterpolatedImage, NewBicubicImageInterpolator())
	InterpolateImage(original, bilinearInterpolatedImage, NewBilinearImageInterpolator())
	InterpolateImage(original, nearestNeighbourInterpolatedImage, NewNearestNeighbourInterpolator())

	writeImage("result/bicubic_interpolated.png", bicubicInterpolatedImage)
	writeImage("result/bilinear_interpolated.png", bilinearInterpolatedImage)
	writeImage("result/nearest_interpolated.png", nearestNeighbourInterpolatedImage)
}

func InterpolateImage(original img.Image, interpolated *img.RGBA, imageInterpolator ImageInterpolator) {
	iw := interpolated.Bounds().Dx()
	ih := interpolated.Bounds().Dy()

	ow := original.Bounds().Dx()
	oh := original.Bounds().Dy()

	dx := float64(ow) / float64(iw)
	dy := float64(oh) / float64(ih)

	for iy := 0; iy < iw; iy++ {
		for ix := 0; ix < iw; ix++ {
			// The float position (fx,fy) in the original image that we want to interpolate
			fx := float64(ix)*dx - 0.5
			fy := float64(iy)*dy - 0.5

			// The actual surrounding four pixels (fx0,fy0) (fx1,fy0) (fx0,fy1) (fx1,fy1) of the pixel to interpolate
			fx0 := int(math.Floor(fx))
			fy0 := int(math.Floor(fy))

			// The percentage distance from (fx0,fy0) to (fx1,fy1)
			fdx := fx - float64(fx0)
			fdy := fy - float64(fy0)

			populateKernel(imageInterpolator.kernelX0, imageInterpolator.kernelY0, imageInterpolator.kernel, fx0, fy0, original)

			color := imageInterpolator.interpolationFunc(fdx, fdy, imageInterpolator.kernel)

			interpolated.Set(ix, iy, color)
		}
	}
}

// https://en.wikipedia.org/wiki/Bicubic_interpolation
func bicubicInterpolateColor(fdx float64, fdy float64, kernel [][]color.Color) color.Color {
	rk := [4][4]float64{}
	gk := [4][4]float64{}
	bk := [4][4]float64{}
	ak := [4][4]float64{}

	// Split RGBA color kernel into kernels for each separate R, G, B, A channel
	var invMax = 1.0 / 0xffff
	for y := 0; y < 4; y++ {
		for x := 0; x < 4; x++ {
			// get RGBA values. The values RGB values has premultiplied alpha.
			// Convert to normal RGB color channels before any interpolation.
			r, g, b, a := kernel[x][y].RGBA()
			alpha := float64(a) * invMax
			invAlpha := 1.0 / alpha
			rk[x][y], gk[x][y], bk[x][y], ak[x][y] = float64(r)*invMax*invAlpha, float64(g)*invMax*invAlpha, float64(b)*invMax*invAlpha, alpha
		}
	}

	r := bicubicInterpolateValue(fdx, fdy, rk)
	g := bicubicInterpolateValue(fdx, fdy, gk)
	b := bicubicInterpolateValue(fdx, fdy, bk)
	a := bicubicInterpolateValue(fdx, fdy, ak)

	// Create interpolated color.
	// The color is supposed to have pre-multiplied alpha so premultiply the alpha to the color.
	conversionConstant := 255.0
	alpha := clampFloat64(a, 0.0, 1.0)
	c := color.RGBA{
		R: uint8(clampFloat64(r*alpha, 0.0, 1.0) * conversionConstant),
		G: uint8(clampFloat64(g*alpha, 0.0, 1.0) * conversionConstant),
		B: uint8(clampFloat64(b*alpha, 0.0, 1.0) * conversionConstant),
		A: uint8(alpha * conversionConstant),
	}

	return c
}

// https://en.wikipedia.org/wiki/Bicubic_interpolation
func bicubicInterpolateValue(fdx float64, fdy float64, kernel [4][4]float64) float64 {
	yKernel := make([]float64, 4)

	for y := 0; y < 4; y++ {
		yKernel[y] = cubicHermiteSplineInterpolateValue(kernel[0][y], kernel[1][y], kernel[2][y], kernel[3][y], fdx)
	}

	v := cubicHermiteSplineInterpolateValue(yKernel[0], yKernel[1], yKernel[2], yKernel[3], fdy)

	return v
}

// https://en.wikipedia.org/wiki/Cubic_Hermite_spline
// https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolation_on_the_unit_interval_with_matched_derivatives_at_endpoints
func cubicHermiteSplineInterpolateValue(pn_1, pn, pn1, pn2, u float64) float64 {
	return 0.5*(((-pn_1+3*pn-3*pn1+pn2)*u+(2*pn_1-5*pn+4*pn1-pn2))*u+(-pn_1+pn1))*u + pn
}

// https://en.wikipedia.org/wiki/Nearest-neighbor_interpolation
func nearestNeighbourInterpolateColor(fdx float64, fdy float64, kernel [][]color.Color) color.Color {
	fdxlow := fdx < 0.5
	fdylow := fdy < 0.5

	c := kernel[0][0]
	if fdxlow && fdylow {
		c = kernel[0][0]
	} else if !fdxlow && fdylow {
		c = kernel[1][0]
	} else if fdxlow && !fdylow {
		c = kernel[0][1]
	} else if !fdxlow && !fdylow {
		c = kernel[1][1]
	}

	return c
}

// https://en.wikipedia.org/wiki/Bilinear_interpolation
func bilinearInterpolateColor(fdx float64, fdy float64, kernel [][]color.Color) color.Color {
	// The known RGBA values for the surrounding four pixels
	f00R, f00G, f00B, f00A := kernel[0][0].RGBA()
	f01R, f01G, f01B, f01A := kernel[0][1].RGBA()
	f10R, f10G, f10B, f10A := kernel[1][0].RGBA()
	f11R, f11G, f11B, f11A := kernel[1][1].RGBA()

	// Bilinear interpolation of (fdx,fdy) for each image channel within the four surrounding pixels
	// with known channel values
	fxyR := bilinearInterpolateValue(fdx, fdy, float64(f00R), float64(f01R), float64(f10R), float64(f11R))
	fxyG := bilinearInterpolateValue(fdx, fdy, float64(f00G), float64(f01G), float64(f10G), float64(f11G))
	fxyB := bilinearInterpolateValue(fdx, fdy, float64(f00B), float64(f01B), float64(f10B), float64(f11B))
	fxyA := bilinearInterpolateValue(fdx, fdy, float64(f00A), float64(f01A), float64(f10A), float64(f11A))

	conversionConstant := 255.0 / float64(0xffff)
	c := color.RGBA{
		R: uint8(fxyR * conversionConstant),
		G: uint8(fxyG * conversionConstant),
		B: uint8(fxyB * conversionConstant),
		A: uint8(fxyA * conversionConstant),
	}

	return c
}

// bilinearInterpolateValue bilineary interpolates a value between the four points (0,0) (1,0) (0,1) (1,1),
// all on the same vertical and horizontal distance from each other (square).
//
// fdx and fdy are the percentage distance in x and y between (0,0) and (1,1) and they are in the range ]0,1].
// f00, f01, f10, f11 are the known values at each point in the square.
//
// f00-------f10
//
//	|  +      |     + is the point (fdx, fdy)
//	|         |       which value we want to interpolate
//	|         |
//
// f01-------f11
//
// https://en.wikipedia.org/wiki/Bilinear_interpolation
func bilinearInterpolateValue(fdx, fdy float64, f00, f01, f10, f11 float64) float64 {
	a00 := f00
	a10 := f10 - f00
	a01 := f01 - f00
	a11 := f11 - f10 - f01 + f00

	fxy := a00 + a10*fdx + a01*fdy + a11*fdx*fdy

	return fxy
}

func populateKernel(kx0 int, ky0 int, kernel [][]color.Color, fx0 int, fy0 int, image img.Image) {
	iw := image.Bounds().Dx()
	ih := image.Bounds().Dy()

	for x := 0; x < len(kernel); x++ {
		for y := 0; y < len(kernel[x]); y++ {
			kernel[x][y] = image.At(clamp(fx0+(x-kx0), 0, iw-1), clamp(fy0+(y-ky0), 0, ih-1))
		}
	}
}

func clamp(x, min, max int) int {
	if x < min {
		return min
	}
	if x > max {
		return max
	}
	return x
}

func clampFloat64(x, min, max float64) float64 {
	if x < min {
		return min
	}
	if x > max {
		return max
	}
	return x
}
