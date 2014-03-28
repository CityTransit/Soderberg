#include <stdio.h>
#include <iostream>
#include <fstream>
#include "frame.h"
#include <stdint.h>
#include <math.h>
#include "gauss.h"

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
        data = new unsigned char[width * height * bitdepth * channels / 8];
        const unsigned int stride = width * bitdepth * channels / 8;

        for(size_t i = 0; i < height; i++) {
            rowPtrs[i] = (png_bytep)data + (height - i - 1)*stride;
        }

        png_read_image(pngPtr, rowPtrs);

        //Contains ALPHA information (we don't want it)
        if(channels == 4) { 
            int new_channels = channels-1;
            unsigned char *new_img = (unsigned char *)calloc(new_channels*width*height, sizeof(char));

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
    unsigned char *new_img = (unsigned char *)calloc(channels*width*height, sizeof(char));
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

bool Frame::applyBilateral(float s, float r)
{
    unsigned char *new_img = (unsigned char *)calloc(channels*width*height, sizeof(char));

    if(!new_img) {
        return false;
    }

    // If they're expecting image values to be [0,1]
    r = (r < 1) ? 255*r : r;
    Gauss* g = new Gauss(r);

    int ix, iy;     //input access
    int kx, ky;     //kernel access
    int endx, endy;
    float pix_gauss, clr_gauss, gauss;

    for (iy=0; iy<height; iy++) {
        for (ix=0; ix<width; ix++) {
            double total[3] = {0};
            double norm[3] = {0};
            int img_pos = iy*width*channels + ix*channels;

            ky = iy - (s-1)/2;
            endx = ix + (s-1)/2;
            endy = iy + (s-1)/2;

            ky = fmax(0, ky);
            endy = fmin(height-1, endy);

            for(; ky<=endy; ky++) {
                kx = ix - (s-1)/2;
                kx = fmax(0, kx); 
                endx = fmin(width-1, endx);
                for(; kx<=endx; kx++) {
                    int knl_pos = ky*width*channels + kx*channels;
                    pix_gauss = gaussian(sqrt(pow((ix - kx), 2) + pow((iy - ky), 2)), s);

                    clr_gauss = g->getGauss(abs(data[img_pos + 0] - data[knl_pos + 0]));
                    gauss = pix_gauss * clr_gauss;
                    total[0] += gauss * data[knl_pos + 0];
                    norm[0] += gauss;

                    clr_gauss = g->getGauss(abs(data[img_pos + 1] - data[knl_pos + 1]));
                    gauss = pix_gauss * clr_gauss;
                    total[1] += gauss * data[knl_pos + 1];
                    norm[1] += gauss;

                    clr_gauss = g->getGauss(abs(data[img_pos + 2] - data[knl_pos + 2]));
                    gauss = pix_gauss * clr_gauss;
                    total[2] += gauss * data[knl_pos + 2];
                    norm[2] += gauss;

                }
            }

            total[0] /= norm[0];
            total[1] /= norm[1];
            total[2] /= norm[2];

            new_img[iy*width*channels + ix*channels + 0] = (unsigned char) fmax(fmin(255, total[0]), 0);
            new_img[iy*width*channels + ix*channels + 1] = (unsigned char) fmax(fmin(255, total[1]), 0);
            new_img[iy*width*channels + ix*channels + 2] = (unsigned char) fmax(fmin(255, total[2]), 0);
        }

    }

    delete data;
    data = new_img;
    
    return true;
}

float Frame::gaussian(float x, float sigma){
    return exp(- x*x / (2*sigma*sigma) ) / (2*M_PI*sigma);
}

bool Frame::applyKernel(Kernel *k)
{
    unsigned char *new_img = (unsigned char *)calloc(channels*width*height, sizeof(char));
    float n = k->get_norm();
    int w = k->get_width();
    int h = k->get_height();
    int offx = (w-1)/2;
    int offy = (h-1)/2;

    if(!new_img) {
        return false;
    }

	int ix, iy;     //input access
    int kx, ky;     //kernel access

	for (iy=0 ; iy<height ; iy++) {
		for (ix=0 ; ix<width ; ix++) {
            int total[3] = {0};
            int pos = iy*width*channels + ix*channels;

            for(ky=0; ky<h; ky++) {

                int imgy = ((iy - offy + ky) + height) % height;
                for(kx=0; kx<w; kx++) {

                    int imgx = (ix - offx + kx + width) % width;
                    //if(iy > 2 && iy < height-3 && ix > 2 && ix < width-3 && kx == 1 && ky == 1 && imgy != iy && imgx != ix) printf("%d %d %d %d\n", imgy, iy, imgx, ix);
                    
                    total[0] += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 0]/n;
                    total[1] += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 1]/n;
                    total[2] += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 2]/n;
                }
            }

            //if(iy == height/2) printf("RGB = %d %d %d\n", total[0], total[1], total[2]);
            new_img[pos + 0] = (unsigned char)fmax(fmin(255, total[0]), 0);
            new_img[pos + 1] = (unsigned char)fmax(fmin(255, total[1]), 0);
            new_img[pos + 2] = (unsigned char)fmax(fmin(255, total[2]), 0);

            //if(new_img[pos] > 255 || new_img[pos+1] > 255 || new_img[pos+2] > 255) printf("RGB = (%d,%d,%d)\n", new_img[iy*width*channels + ix*channels + 0], new_img[iy*width*channels + ix*channels + 0], new_img[iy*width*channels + ix*channels + 0]);
		}

	}
    delete data;
    data = new_img;
    
    return true;
}

bool Frame::applyTwoToneKernel(Kernel *k)
{
    unsigned char *new_img = (unsigned char *)calloc(channels*width*height, sizeof(char));

    int w = k->get_width();
    int h = k->get_height();
    int offx = (w-1)/2;
    int offy = (h-1)/2;

    if(!new_img) {
        return false;
    }

    int ix, iy;     //input access
    int kx, ky;     //kernel access

    for (iy=0 ; iy<height ; iy++) {
        for (ix=0 ; ix<width ; ix++) {
            float total = 0;
            int pos = iy*width*channels + ix*channels;

            for(ky=0; ky<h; ky++) {
                int imgy = ((iy - offy + ky) + height) % height;
                for(kx=0; kx<w; kx++) {
                    int imgx = (ix - offx + kx + width) % width;

                    float tmp = 0;
                    tmp += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 0];
                    tmp += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 1];
                    tmp += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 2];
                    tmp/=3;
                    total+=tmp;
                }
            }

            if(total <= 0) total = 255;
            else total = 0;

            new_img[pos + 0] = total;
            new_img[pos + 1] = total;
            new_img[pos + 2] = total;

        }

    }
    delete data;
    data = new_img;
    
    return true;
}


bool Frame::applyDoG(float sigma, float k)
{
    Kernel *k1 = Kernel::generateGaussian(sigma);
    Kernel *k2 = Kernel::generateGaussian(sigma * k);
    unsigned char *new_img = (unsigned char *)calloc(channels*width*height, sizeof(char));
    float n1 = k1->get_norm();
    int w1 = k1->get_width();
    int h1 = k1->get_height();
    int offx1 = (w1-1)/2;
    int offy1 = (h1-1)/2;
    

    float n2 = k2->get_norm();
    int w2 = k2->get_width();
    int h2 = k2->get_height();
    int offx2 = (w2-1)/2;
    int offy2 = (h2-1)/2;

    if(w2 < w1 || h2 < h1){
        perror("Bad kernel sizes applyDoG()");
        return false;
    }

    if(!new_img) {
        return false;
    }

    int ix, iy;     //input access
    int kx, ky;     //kernel access

    for (iy=0 ; iy<height ; iy++) {
        for (ix=0 ; ix<width ; ix++) {
            float total1 = 0;
            float total2 = 0;
            int pos = iy*width*channels + ix*channels;

            for(ky=0; ky<h2; ky++) {
                int imgy1 = ((iy - offy1 + ky) + height) % height;
                int imgy2 = ((iy - offy2 + ky) + height) % height;
                for(kx=0; kx<w2; kx++) {
                    if(ky<h1 && kx<w1){
                        int imgx1 = (ix - offx1 + kx + width) % width;
                        int pos2 = (imgy1*width*channels + imgx1*channels);

                        float tmp = data[pos2 + 0] + data[pos2 + 1] + data[pos2 + 2];
                        tmp *= k1->get(ky*w1 + kx);
                        tmp/=3 * n1;

                        total1+=tmp;
                    }
                    int imgx2 = (ix - offx2 + kx + width) % width;
                    int pos2 = (imgy2*width*channels + imgx2*channels);

                    float tmp = data[pos2 + 0] + data[pos2 + 1] + data[pos2 + 2];
                    tmp *= k2->get(ky*w2 + kx);
                    tmp/=3 * n2;

                    total2+=tmp;
                }
            }
            
            float dog_val;// = threshold(total1, total2, 0);

            float p = 0.2;
            p *= 255;
            dog_val = (1 + p) * total1 - p * total2;
            
            dog_val = (unsigned char)fmax(fmin(255, dog_val), 0);

            new_img[pos + 0] =
            new_img[pos + 1] =
            new_img[pos + 2] = dog_val;

        }

    }
    delete data;
    data = new_img;
    
    return true;
}

float Frame::threshold(float a, float b, float t){
    float res;
    if((a -b) > t){
        res = 255 - (a - b);
    } else {
        res = a - b;
    }
    return res;
}
