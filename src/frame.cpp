#include <stdio.h>
#include <iostream>
#include <fstream>
#include "frame.h"
#include <stdint.h>
#include <math.h>

void Frame::userReadData(png_structp pngPtr, png_bytep data, png_size_t length)
{
    png_voidp a = png_get_io_ptr(pngPtr);

    ((std::istream *)a)->read((char *)data, length);
}

bool Frame::validate(std::istream &src)
{
    png_byte sig[PNGSIGSIZE];

    src.read((char *)sig, PNGSIGSIZE);

    if(!src.good()) {
        return false;
    }

    return (png_sig_cmp(sig, 0, PNGSIGSIZE) == 0);
}

bool Frame::open(const char *fname) 
{
    std::filebuf img;

    if(img.open(fname, std::ios::in)) {
        std::istream smg(&img);

        if(!validate(smg)) {
            fprintf(stderr, "ERROR: %s is not a valid png file.\n", fname);
            return false;
        }

        pngPtr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!pngPtr) {
            fprintf(stderr, "ERROR: Unable to initialize png struct.\n");
            return false;
        }

        infoPtr = png_create_info_struct(pngPtr);
        if(!infoPtr) {
            fprintf(stderr, "ERROR: Unable to initialize ino struct.\n");
            png_destroy_read_struct(&pngPtr, (png_infopp)0, (png_infopp)0);
            return false;
        }

        rowPtrs = NULL;
        data = NULL;

        if(setjmp(png_jmpbuf(pngPtr))) {
            png_destroy_read_struct(&pngPtr, &infoPtr, (png_infopp)0);
            if(rowPtrs) {
                delete [] rowPtrs;
            }
            if(data) {
                delete [] data;
            }

            fprintf(stderr, "ERROR: An error occurred while reading %s.\n", fname);

            return false;
        }

        png_set_read_fn(pngPtr, (png_voidp)&smg, userReadData);
        png_set_sig_bytes(pngPtr, PNGSIGSIZE);
        png_read_info(pngPtr, infoPtr);

        width = png_get_image_width(pngPtr, infoPtr);
        height = png_get_image_height(pngPtr, infoPtr);

        bitdepth = png_get_bit_depth(pngPtr, infoPtr);
        channels = png_get_channels(pngPtr, infoPtr);
        colourtype = png_get_color_type(pngPtr, infoPtr);

        switch(colourtype) {
            case PNG_COLOR_TYPE_PALETTE:
                png_set_palette_to_rgb(pngPtr);
                channels = 3;
                break;
            case PNG_COLOR_TYPE_GRAY:
                if(bitdepth < 8) {
                    png_set_expand_gray_1_2_4_to_8(pngPtr);
                }
                bitdepth = 8;
                break;
        }

        if(png_get_valid(pngPtr, infoPtr, PNG_INFO_tRNS)) {
            png_set_tRNS_to_alpha(pngPtr);
            channels += 1;
        }

        if(bitdepth == 16) {
            png_set_strip_16(pngPtr);
        }

        rowPtrs = new png_bytep[height];
        data = new char[width * height * bitdepth * channels / 8];
        const unsigned int stride = width * bitdepth * channels / 8;

        for(size_t i = 0; i < height; i++) {
            rowPtrs[i] = (png_bytep)data + (height - i - 1)*stride;
        }

        png_read_image(pngPtr, rowPtrs);

        //Contains ALPHA information (we don't want it)
        if(channels == 4) { 
            int new_channels = channels-1;
            char *new_img = (char *)calloc(new_channels*width*height, sizeof(char));

            if(!new_img) {
                fprintf(stderr, "WARNING: Image still contains ALPHA information.\n");
                return false;
            }
            
            // Write image data
            int x, y;
            for (y=0 ; y<height ; y++) {
                for (x=0 ; x<width ; x++) {
                    new_img[y*width*new_channels + x*new_channels + 0] = data[y*width*channels + x*channels + 0];
                    new_img[y*width*new_channels + x*new_channels + 1] = data[y*width*channels + x*channels + 1];
                    new_img[y*width*new_channels + x*new_channels + 2] = data[y*width*channels + x*channels + 2];
                }
            }

            delete data;
            data = new_img;
            channels = new_channels;
        }
            
        //delete [] (png_bytep)rowPtrs;
        //png_destroy_read_struct(&pngPtr, &infoPtr, (png_infopp)0);
    }
    else {
        fprintf(stderr, "ERROR: Unable to open file %s.\n", fname);
    }

    return true;
}

// Deconstructor
Frame::~Frame()
{
    delete [] (png_bytep)rowPtrs;
    png_destroy_read_struct(&pngPtr, &infoPtr, (png_infopp)0);
    //delete image;
}

bool Frame::save(const char *filename, char *title)
{
	int code = 0;
	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;
	
	// Open file for writing (binary mode)
	fp = fopen(filename, "wb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s for writing\n", filename);
		code = 1;
		goto finalise;
	}
    
	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fprintf(stderr, "ERROR: Could not allocate write struct\n");
		code = 1;
		goto finalise;
	}
    
	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fprintf(stderr, "ERROR: Could not allocate info struct\n");
		code = 1;
		goto finalise;
	}
    
	// Setup Exception handling
	if (setjmp(png_jmpbuf(png_ptr))) {
		fprintf(stderr, "ERROR: Error during png creation\n");
		code = 1;
		goto finalise;
	}
    
	png_init_io(png_ptr, fp);
    
	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    
	// Set title
	if (title != NULL) {
		png_text title_text;
		title_text.compression = PNG_TEXT_COMPRESSION_NONE;
		title_text.key = "Title";
		title_text.text = title;
		png_set_text(png_ptr, info_ptr, &title_text, 1);
	}
    
	png_write_info(png_ptr, info_ptr);
    
	// Allocate memory for one row (4 bytes per pixel - RGBA)
	row = (png_bytep) malloc(3 * width * sizeof(png_byte));
    
	// Write image data
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++) {
			row[x*3 + 0] = data[(height-y-1)*width*channels + x*channels + 0];
			row[x*3 + 1] = data[(height-y-1)*width*channels + x*channels + 1];
			row[x*3 + 2] = data[(height-y-1)*width*channels + x*channels + 2];
		}
		png_write_row(png_ptr, row);
	}
    
	// End write
	png_write_end(png_ptr, NULL);
    
finalise:
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);
    
	return code;
}

bool Frame::flip()
{
    char *new_img = (char *)calloc(channels*width*height, sizeof(char));
	//row = (png_bytep) malloc(3 * width * sizeof(png_byte));
    if(!new_img) {
        return false;
    }
    
	// Write image data
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++) {
			new_img[y*width*channels + x*channels + 0] = data[(height-y-1)*width*channels + x*channels + 0];
			new_img[y*width*channels + x*channels + 1] = data[(height-y-1)*width*channels + x*channels + 1];
			new_img[y*width*channels + x*channels + 2] = data[(height-y-1)*width*channels + x*channels + 2];
		}
	}

    delete data;
    data = new_img;
    
    return true;
}
