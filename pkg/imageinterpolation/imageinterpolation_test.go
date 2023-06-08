package imageinterpolation

import (
	"fmt"
	img "image"
	"image/png"
	"os"
	"path/filepath"
	"strconv"
	"testing"
)

func Test_BicubicUpscale(t *testing.T) {
	original, err := readImage("../../resources/wiki.png")
	if err != nil {
		panic(err)
	}

	upscale := 50
	newWidth := original.Bounds().Dx() * upscale
	newHeight := original.Bounds().Dy() * upscale
	bicubicInterpolatedImage := img.NewRGBA(img.Rect(0, 0, newWidth, newHeight))

	InterpolateImage(original, bicubicInterpolatedImage, NewBicubicImageInterpolator())

	writeImage("../../result/bicubic_interpolated.png", bicubicInterpolatedImage)
}

func Test_BilinearUpscale(t *testing.T) {
	original, err := readImage("../../resources/wiki.png")
	if err != nil {
		panic(err)
	}

	upscale := 50
	newWidth := original.Bounds().Dx() * upscale
	newHeight := original.Bounds().Dy() * upscale
	bilinearInterpolatedImage := img.NewRGBA(img.Rect(0, 0, newWidth, newHeight))

	InterpolateImage(original, bilinearInterpolatedImage, NewBilinearImageInterpolator())

	writeImage("../../result/bilinear_interpolated.png", bilinearInterpolatedImage)
}

func Test_NearestNeighbourUpscale(t *testing.T) {
	original, err := readImage("../../resources/wiki.png")
	if err != nil {
		panic(err)
	}

	upscale := 50
	newWidth := original.Bounds().Dx() * upscale
	newHeight := original.Bounds().Dy() * upscale

	nearestNeighbourInterpolatedImage := img.NewRGBA(img.Rect(0, 0, newWidth, newHeight))

	InterpolateImage(original, nearestNeighbourInterpolatedImage, NewNearestNeighbourInterpolator())
	writeImage("../../result/nearest_interpolated.png", nearestNeighbourInterpolatedImage)
}

func Test_Resize(t *testing.T) {
	boxOriginal, err := readImage("../../resources/library.png") // original image is 800x800 pixels
	if err != nil {
		panic(err)
	}

	resizedImage0 := ResizeImage(boxOriginal, boxOriginal.Bounds().Dx(), boxOriginal.Bounds().Dy(), NewBicubicImageInterpolator())
	resizedImage1 := ResizeImage(boxOriginal, boxOriginal.Bounds().Dx()/4, boxOriginal.Bounds().Dy()/4, NewBicubicImageInterpolator())
	resizedImage2 := ResizeImage(boxOriginal, 473, 0, NewBicubicImageInterpolator())
	resizedImage3 := ResizeImage(boxOriginal, 0, 473, NewBicubicImageInterpolator())
	resizedImage4 := ResizeImage(boxOriginal, 100, 200, NewBicubicImageInterpolator())
	resizedImage5 := ResizeImage(boxOriginal, 200, 100, NewBicubicImageInterpolator())
	resizedImage6 := ResizeImage(boxOriginal, 473, 967, NewBicubicImageInterpolator())
	resizedImage7 := ResizeImage(boxOriginal, 967, 473, NewBicubicImageInterpolator())
	resizedImage8 := ResizeImage(boxOriginal, 1367, 813, NewBicubicImageInterpolator())

	writeImage("../../result/resized_0.png", resizedImage0)
	writeImage("../../result/resized_1.png", resizedImage1)
	writeImage("../../result/resized_2.png", resizedImage2)
	writeImage("../../result/resized_3.png", resizedImage3)
	writeImage("../../result/resized_4.png", resizedImage4)
	writeImage("../../result/resized_5.png", resizedImage5)
	writeImage("../../result/resized_6.png", resizedImage6)
	writeImage("../../result/resized_7.png", resizedImage7)
	writeImage("../../result/resized_8.png", resizedImage8)
}

func Test_BoxInterpolation(t *testing.T) {
	boxOriginal, err := readImage("../../resources/library.png") // original image is 800x800 pixels
	if err != nil {
		panic(err)
	}

	boxInterpolatedImage := boxOriginal
	for boxIteration := 0; boxIteration < 4; boxIteration++ {
		boxInterpolatedImage = boxInterpolation(2, 2, true, boxInterpolatedImage)
		writeImage("../../result/box_interpolated_"+strconv.Itoa(boxIteration+1)+".png", boxInterpolatedImage)
	}
}

func readImage(filePath string) (img.Image, error) {
	f, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	image, _, err := img.Decode(f)
	return image, err
}

func writeImage(filename string, image img.Image) {
	parentPath := filepath.Dir(filename)
	os.MkdirAll(parentPath, os.ModePerm)

	f, err := os.Create(filename)
	if err != nil {
		fmt.Println("Oups, no files for you today.")
		os.Exit(1)
	}
	defer f.Close()

	// Encode to `PNG` with `DefaultCompression` level then save to file
	err = png.Encode(f, image)
	if err != nil {
		fmt.Println("Oups, no image encode for you today.")
		os.Exit(1)
	}
}
