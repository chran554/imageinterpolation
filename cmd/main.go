package main

import (
	"fmt"
	"image/color"
	"math"
	"os"
	img "image"
	_ "image/jpeg"
	"image/png"
)

func main() {
	fmt.Println("Up and running...")

	original, err := getImageFromFilePath("resources/original.png")
	if err != nil {
		panic(err)
	}

	width := original.Bounds().Dx() * 10
	height := original.Bounds().Dy() * 10

	bilinearInterpolatedImage := img.NewRGBA(img.Rect(0, 0, width, height))


	bilinearInterpolation(original, bilinearInterpolatedImage)
}

func bilinearInterpolation(original img.Image, interpolated *img.RGBA) {
	iw := interpolated.Bounds().Dx()
	ih := interpolated.Bounds().Dy()

	dx := float64(original.Bounds().Dx()) / float64(iw)
	dy := float64(original.Bounds().Dy()) / float64(ih)

	for iy := 0; iy < iw; iy++ {
		for ix := 0; ix < iw; ix++ {
			fx := float64(ix) * dx
			fy := float64(iy) * dy

			fx0 := math.Floor(fx)
			fy0 := math.Floor(fy)
			fx1 := fx0 +1
			fy1 := fy0 +1

			fdx := fx -fx0
			fdy := fy -fy0

			f00 := original.At(int(fx0), int(fy0))
			f01 := original.At(int(fx0), int(fy1))
			f10 := original.At(int(fx1), int(fy0))
			f11 := original.At(int(fx1), int(fy1))

			f00R, f00G, f00B, f00A := f00.RGBA()
			f01R, f01G, f01B, f01A := f01.RGBA()
			f10R, f10G, f10B, f10A := f10.RGBA()
			f11R, f11G, f11B, f11A := f11.RGBA()

			fxyR := bilinearInterpolate(fdx, fdy, float64(f00R), float64(f01R), float64(f10R), float64(f11R))
			fxyG := bilinearInterpolate(fdx, fdy, float64(f00G), float64(f01G), float64(f10G), float64(f11G))
			fxyB := bilinearInterpolate(fdx, fdy, float64(f00B), float64(f01B), float64(f10B), float64(f11B))
			fxyA := bilinearInterpolate(fdx, fdy, float64(f00A), float64(f01A), float64(f10A), float64(f11A))

			interpolated.Set(ix, iy, color.Color)
		}
	}
}

func bilinearInterpolate(fdx, fdy float64, f00, f01, f10, f11 float64) float64 {
	a00 := f00
	a10 := f10 - f00
	a01 := f01 - f00
	a11 := f11 - f10 - f01 + f00

	fxy := a00 + a10 * fdx + a01 * fdy + a11 *fdx * fdy

	return fxy
}

func getImageFromFilePath(filePath string) (img.Image, error) {
	f, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	image, _, err := img.Decode(f)
	return image, err
}

