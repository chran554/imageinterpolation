package main

import (
	"fmt"
	img "image"
	"image/png"
	"os"
	"path/filepath"
)

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
