package imageinterpolation

import (
	img "image"
	"image/color"
	_ "image/jpeg"
	_ "image/png"
	"math"
)

type ImageInterpolator struct {
	kernelX0          int
	kernelY0          int
	kernel            [][]*linearRGBA
	interpolationFunc func(float64, float64, [][]*linearRGBA) color.Color
}

func NewBicubicImageInterpolator() ImageInterpolator {
	return ImageInterpolator{
		kernelX0:          1,
		kernelY0:          1,
		kernel:            createColorKernel(4, 4),
		interpolationFunc: bicubicInterpolateColor,
	}
}

func NewBilinearImageInterpolator() ImageInterpolator {
	return ImageInterpolator{
		kernelX0:          0,
		kernelY0:          0,
		kernel:            createColorKernel(2, 2),
		interpolationFunc: bilinearInterpolateColor,
	}
}

func NewNearestNeighbourInterpolator() ImageInterpolator {
	return ImageInterpolator{
		kernelX0:          0,
		kernelY0:          0,
		kernel:            createColorKernel(2, 2),
		interpolationFunc: nearestNeighbourInterpolateColor,
	}
}

// ResizeImage resizes an image to new dimensions (width and height).
//
// One of the parameters newWidth and newHeight can be of value 0, but not both.
// Value 0 implies that the resize should calculate that dimension of the resized image to keep
// the original image aspect ratio between width and height.
func ResizeImage(original img.Image, newWidth, newHeight int, imageInterpolator ImageInterpolator) *img.RGBA {
	var resizedImage *img.RGBA

	if (newWidth == 0) && (newHeight == 0) {
		newWidth = original.Bounds().Dx()
		newHeight = original.Bounds().Dy()
	} else if newWidth == 0 {
		newWidth = int(float64(original.Bounds().Dx()) * (float64(newHeight) / float64(original.Bounds().Dy())))
	} else if newHeight == 0 {
		newHeight = int(float64(original.Bounds().Dy()) * (float64(newWidth) / float64(original.Bounds().Dx())))
	}

	resizedImage = img.NewRGBA(img.Rect(0, 0, newWidth, newHeight))

	// If we have matching bounds on our original image to the desired dimensions then return a copy of the original image.
	if original.Bounds().Dx() == newWidth && original.Bounds().Dy() == newHeight {
		for y := 0; y < newHeight; y++ {
			for x := 0; x < newWidth; x++ {
				resizedImage.Set(x, y, original.At(x, y))
			}
		}
		return resizedImage
	}

	// Scale down original image to at least size of resized image (both dimensions)
	// and at most double the size (both dimensions).
	// If any original dimension is less than new size then it is ok but that dimension is left untouched.
	downscale := true
	for downscale {
		ow := original.Bounds().Dx()
		oh := original.Bounds().Dy()

		scaleWidth := 1
		if (ow / newWidth) >= 2 {
			scaleWidth = 2
		}

		scaleHeight := 1
		if (oh / newHeight) >= 2 {
			scaleHeight = 2
		}

		downscale = (scaleWidth > 1) || (scaleHeight > 1)
		if downscale {
			original = boxInterpolation(scaleWidth, scaleHeight, true, original)
		}
	}

	// If we have matching bounds on our original image to the desired dimensions then return a copy of the original image.
	if original.Bounds().Dx() == newWidth && original.Bounds().Dy() == newHeight {
		for y := 0; y < newHeight; y++ {
			for x := 0; x < newWidth; x++ {
				resizedImage.Set(x, y, original.At(x, y))
			}
		}
		return resizedImage
	}

	// Interpolate the (perhaps downscaled) original to the resized image.
	InterpolateImage(original, resizedImage, imageInterpolator)

	return resizedImage
}

func InterpolateImage(original img.Image, interpolated *img.RGBA, imageInterpolator ImageInterpolator) {
	iw := interpolated.Bounds().Dx()
	ih := interpolated.Bounds().Dy()

	ow := original.Bounds().Dx()
	oh := original.Bounds().Dy()

	dx := float64(ow) / float64(iw)
	dy := float64(oh) / float64(ih)

	for iy := 0; iy < ih; iy++ {
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

			populateKernel(imageInterpolator.kernelX0, imageInterpolator.kernelY0, imageInterpolator.kernel, fx0, fy0, original, 2.2)

			color := imageInterpolator.interpolationFunc(fdx, fdy, imageInterpolator.kernel)

			interpolated.Set(ix, iy, color)
		}
	}
}

// boxInterpolation reduces an image by a horizontal and vertical factor.
// Destination image has a reduced width of (1.0/horizontalScaleDivisor)*original.width and a reduced height of (1.0/verticalScaleDivisor)*original.height.
// Any scale factor of 1 will leave that image dimension unaffected.
// Any scale factor of 0 or less is prohibited.
//
// Interpolation is done through taking a box of pixels (of size horizontalScaleDivisor x verticalScaleDivisor) in the original image,
// calculate an average color and put that color on a pixel in the destination image.
//
// If original image dimensions width and height is not perfectly divisible by the factors
// any remainder pixels at the right or bottom will be discarded unless you choose to set parameter perfect to true.
// Perfect will include remainder pixels in destination image.
func boxInterpolation(horizontalScaleDivisor, verticalScaleDivisor int, perfect bool, original img.Image) (interpolated *img.RGBA) {
	kernel := createColorKernel(horizontalScaleDivisor, verticalScaleDivisor)

	ow := original.Bounds().Dx()
	oh := original.Bounds().Dy()

	iw := ow / horizontalScaleDivisor
	ih := oh / verticalScaleDivisor

	// Prepare original image for boxing.
	horizontalRemainder := math.Remainder(float64(ow), float64(horizontalScaleDivisor))
	verticalRemainder := math.Remainder(float64(oh), float64(verticalScaleDivisor))
	if perfect && ((horizontalRemainder > 0.0) || (verticalRemainder > 0.0)) {
		interpolated := img.NewRGBA(img.Rect(0, 0, ow-int(horizontalRemainder), oh-int(verticalRemainder)))
		InterpolateImage(original, interpolated, NewBicubicImageInterpolator())
		original = interpolated
	}

	interpolated = img.NewRGBA(img.Rect(0, 0, iw, ih))

	for iy := 0; iy < ih; iy++ {
		for ix := 0; ix < iw; ix++ {
			populateKernel(0, 0, kernel, ix*horizontalScaleDivisor, iy*verticalScaleDivisor, original, 2.2)
			c := boxInterpolate(kernel)
			interpolated.Set(ix, iy, c)
		}
	}

	return interpolated
}

func boxInterpolate(kernel [][]*linearRGBA) color.Color {
	var r, g, b, a uint32 = 0, 0, 0, 0
	for x := 0; x < len(kernel); x++ {
		for y := 0; y < len(kernel[x]); y++ {
			// TODO de-premultiply alpha
			dr, dg, db, da := kernel[x][y].RGBA()
			r += dr
			g += dg
			b += db
			a += da
		}
	}

	conversion := 255.0 / float64(len(kernel)*len(kernel[0])*0xffff)
	// TODO premultiply alpha back again
	return color.RGBA{R: uint8(float64(r) * conversion), G: uint8(float64(g) * conversion), B: uint8(float64(b) * conversion), A: uint8(float64(a) * conversion)}
}

// https://en.wikipedia.org/wiki/Bicubic_interpolation
func bicubicInterpolateColor(fdx float64, fdy float64, kernel [][]*linearRGBA) color.Color {
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
func nearestNeighbourInterpolateColor(fdx float64, fdy float64, kernel [][]*linearRGBA) color.Color {
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
func bilinearInterpolateColor(fdx float64, fdy float64, kernel [][]*linearRGBA) color.Color {
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

func createColorKernel(width int, height int) [][]*linearRGBA {
	kernel := make([][]*linearRGBA, width)

	for columnIndex := range kernel {
		kernel[columnIndex] = make([]*linearRGBA, height)
	}

	return kernel
}

func populateKernel(kx0 int, ky0 int, kernel [][]*linearRGBA, fx0 int, fy0 int, image img.Image, gamma float64) {
	iw := image.Bounds().Dx()
	ih := image.Bounds().Dy()

	for x := 0; x < len(kernel); x++ {
		for y := 0; y < len(kernel[x]); y++ {
			pixelColor := image.At(clamp(fx0+(x-kx0), 0, iw-1), clamp(fy0+(y-ky0), 0, ih-1))
			kernel[x][y] = newLinearRGBA(pixelColor, gamma)
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

// linearRGBA
//
// https://blog.johnnovak.net/2016/09/21/what-every-coder-should-know-about-gamma/
type linearRGBA struct {
	gamma      float64
	c          *color.Color
	R, G, B, A float64
}

// RGBA
//
// https://blog.johnnovak.net/2016/09/21/what-every-coder-should-know-about-gamma/
func (lRGBA *linearRGBA) RGBA() (r, g, b, a uint32) {
	if lRGBA.c == nil || true { // TODO check if c is already calculated (not nil), if so, do not redo this calculation
		invGamma := 1.0 / lRGBA.gamma
		premultiplyAlphaConstant := lRGBA.A
		c := color.Color(color.RGBA{
			R: uint8(clampFloat64(gammaCalculation(lRGBA.R*premultiplyAlphaConstant, invGamma)*255.0, 0.0, 255.0)),
			G: uint8(clampFloat64(gammaCalculation(lRGBA.G*premultiplyAlphaConstant, invGamma)*255.0, 0.0, 255.0)),
			B: uint8(clampFloat64(gammaCalculation(lRGBA.B*premultiplyAlphaConstant, invGamma)*255.0, 0.0, 255.0)),
			A: uint8(clampFloat64(lRGBA.A*255.0, 0.0, 255.0)),
		})

		lRGBA.c = &c
	}

	return (*lRGBA.c).RGBA()
}

// newLinearRGBA
//
// https://blog.johnnovak.net/2016/09/21/what-every-coder-should-know-about-gamma/
func newLinearRGBA(c color.Color, gamma float64) *linearRGBA {
	r, g, b, a := c.RGBA()

	const conversionConstant = 1.0 / float64(0xffff)
	depremultiplyAlphaConstant := 1.0 / (float64(a) * conversionConstant)
	linearColor := &linearRGBA{
		gamma: gamma,
		c:     &c,
		R:     gammaCalculation(float64(r)*conversionConstant*depremultiplyAlphaConstant, gamma),
		G:     gammaCalculation(float64(g)*conversionConstant*depremultiplyAlphaConstant, gamma),
		B:     gammaCalculation(float64(b)*conversionConstant*depremultiplyAlphaConstant, gamma),
		A:     float64(a) * conversionConstant,
	}

	return linearColor
}

func gammaCalculation(value float64, gamma float64) float64 {
	return math.Pow(value, gamma)
}
