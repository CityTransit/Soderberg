#include <stdio.h>
#include <iostream>
#include <fstream>
#include "frame.h"
#include <stdint.h>
#include <math.h>
#include "gauss.h"

#define LENGTH_2(v) sqrt(v[0]*v[0] + v[1]*v[1])
#define LENGTH_3(v) sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
#define DOT_2(u,v) (u[0]*v[0] + u[1]*v[1])
#define DOT_3(u,v) (u[0]*v[0] + u[1]*v[1] + u[2]*v[2])
#define MAT_MULT_22(u,v) {u[0]*v[0] + u[1]*v[2], u[0]*v[1] + u[1]*v[3], u[2]*v[0] + u[3]*v[2], u[2]*v[1] + u[3]*v[3]}
#define MAT_MULT_21(u,v) {u[0]*v[0] + u[1]*v[1], u[2]*v[0] + u[3]*v[1]}

Frame::Frame(Frame *copy)
{
    this->width = copy->width;
    this->height = copy->height;
    this->bitdepth = copy->bitdepth;
    this->channels = copy->channels;
    this->colourtype = copy->colourtype;

    memcpy(this->data, copy->data, (channels*width*height*sizeof(char)));
}

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
    //delete [] (png_bytep)rowPtrs;
    //png_destroy_read_struct(&pngPtr, &infoPtr, (png_infopp)0);
    //if(data) free(data);
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

/* flip()
 *
 * PURPOSE
 * Flips the image
 * 
 * INPUT:
 * n/a
 * 
 * EXPECTATIONS:
 * n/a
 * 
 * NOTE: 
 * n/a
 * 
 * OUTPUT:
 * Applies filter to global image: char* data 
 * 
 */

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

/* applyBilateral(float S, float R)
 * 
 * PURPOSE:
 * Applies a bilateral filter to the image
 * 
 * INPUT:
 * S represents the spatial distance of the Gaussian kernel
 * R represents the standard deviation of pixel values
 * 
 * EXPECTATIONS:
 * S should be a positive natural number, with greater values giving increased smoothness
 * R should be in range [0...1], with lower values giving increased blur
 *
 * NOTE: 
 * Larger S values have a negative impact on performance
 * 
 * OUTPUT:
 * Applies filter to global image: char* data 
 * 
 */

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

    //Iterate over the original image
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

            //Iterate over the Kernel
            for(; ky<=endy; ky++) {
                kx = ix - (s-1)/2;
                kx = fmax(0, kx); 
                endx = fmin(width-1, endx);
                for(; kx<=endx; kx++) {

                    //Apply the bilateral
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

/*
 * Applys a temporal bilateral filter to the current frame.
 *
 * s - Spatial gaussian standard deviation
 * r - Pixel value gaussian standard deviation
 * t - Temporal gaussian standard deviation
 * *f - Pointer to an array of Frames
 * fn - Length of the pointer or arrays
 *
 */

bool Frame::applyTemporalBilateral(float s, float r, float t, Frame *f, int fn)
{
    unsigned char *new_img = (unsigned char *)calloc(channels*width*height, sizeof(char));

    if(!new_img) {
        return false;
    }

    // If they're expecting image values to be [0,1]
    r = (r < 1) ? 255*r : r;

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

            // Apply spatial and pixel value filters
            for(; ky<=endy; ky++) {
                kx = ix - (s-1)/2;
                kx = fmax(0, kx); 
                endx = fmin(width-1, endx);
                for(; kx<=endx; kx++) {
                    int knl_pos = ky*width*channels + kx*channels;
                    pix_gauss = gaussian(sqrt(pow((ix - kx), 2) + pow((iy - ky), 2)), s);

                    clr_gauss = gaussian(abs(data[img_pos + 0] - data[knl_pos + 0]), r);
                    gauss = pix_gauss * clr_gauss;
                    total[0] += gauss * data[knl_pos + 0];
                    norm[0] += gauss;

                    clr_gauss = gaussian(abs(data[img_pos + 1] - data[knl_pos + 1]), r);
                    gauss = pix_gauss * clr_gauss;
                    total[1] += gauss * data[knl_pos + 1];
                    norm[1] += gauss;

                    clr_gauss = gaussian(abs(data[img_pos + 2] - data[knl_pos + 2]), r);
                    gauss = pix_gauss * clr_gauss;
                    total[2] += gauss * data[knl_pos + 2];
                    norm[2] += gauss;

                }
            }

            // Apply temporal filter
            for(int p=0; p<fn; p++) {
                gauss = gaussian( abs(data[img_pos + 0] - f[p].get(ix, iy, 0)), t);
                total[0] += gauss * f[p].get(ix, iy, 0);
                norm[0] += gauss;

                gauss = gaussian( abs(data[img_pos + 1] - f[p].get(ix, iy, 1)), t);
                total[1] += gauss * f[p].get(ix, iy, 1);
                norm[1] += gauss;

                gauss = gaussian( abs(data[img_pos + 2] - f[p].get(ix, iy, 2)), t);
                total[2] += gauss * f[p].get(ix, iy, 2);
                norm[2] += gauss;
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

/* applyKernel(Kernel *k)
 * 
 * PURPOSE:
 * Applies any pre-generated Kernel to the image
 * 
 * INPUT:
 * A pre-generated kernel of type Kernel*
 * 
 * EXPECTATIONS:
 * k should have equal height and width
  *
 * NOTE: 
 * Larger kernels will have a negative impact on performance
 * 
 * OUTPUT:
 * Applies kernel to global image: char* data 
 * 
 */

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
                    
                    total[0] += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 0]/n;
                    total[1] += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 1]/n;
                    total[2] += k->get(ky*w + kx) * data[(imgy*width*channels + imgx*channels) + 2]/n;
                }
            }

            new_img[pos + 0] = (unsigned char)fmax(fmin(255, total[0]), 0);
            new_img[pos + 1] = (unsigned char)fmax(fmin(255, total[1]), 0);
            new_img[pos + 2] = (unsigned char)fmax(fmin(255, total[2]), 0);
		}
	}
    delete data;
    data = new_img;
    
    return true;
}

/* applyTwoToneKernel(Kernel *k)
 * 
 * PURPOSE:
 * Applies any pre-generated Kernel to the image, and thresholds the result to black/white
 * 
 * INPUT:
 * A pre-generated kernel of type Kernel*
 * 
 * EXPECTATIONS:
 * k should have equal height and width
 *
 * NOTE: 
 * Larger kernels will have a negative impact on performance
 * 
 * OUTPUT:
 * Applies kernel to global image: char* data 
 * 
 */

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

/* applyDoG(float sigma, float k)
 * 
 * PURPOSE:
 * Applies a DoG filter to the image, thresholds the result in B&W
 * 
 * INPUT:
 * float sigma: The size of Kernel A
 * float k: The scalar with which to augment the size of Kernel B
 * 
 * EXPECTATIONS:
 * k should be 1 or above, (1 will have no effect)
 * It is preferable to have the image pre-filtered with a bilateral filter
 *
 * NOTE: 
 * Larger sigma or k values will have a negative impact on performance
 * 
 * OUTPUT:
 * Applies DoG edge finding to global image: char* data 
 * 
 */

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
            //Threshold DoG
            float dog_val = threshold(total1, total2, 0);
            dog_val = (unsigned char)fmax(fmin(255, dog_val), 0);

            new_img[pos + 0] = 
                new_img[pos + 1] = 
                new_img[pos + 2] = dog_val;

        }

    }
    delete data;
    delete k1;
    delete k2;
    data = new_img;
    
    return true;
}

/* applyXDoG(float sigma, float k)
 * 
 * PURPOSE:
 * Applies a XDoG filter to the image
 * XDoG is a DoG filter augmented with values for customization
 * and smoother color transitions.
 * Thresholds the result in B&W
 * 
 * INPUT:
 * float sigma: The size of Kernel A
 * float k: The scalar with which to augment the size of Kernel B
 * 
 * EXPECTATIONS:
 * k should be 1 or above, (1 will have no effect)
 * It is preferable to have the image pre-filtered with a bilateral filter
 *
 * NOTE: 
 * Larger sigma or k values will have a negative impact on performance
 * 
 * OUTPUT:
 * Applies XDoG edge finding to global image: char* data 
 * 
 */

bool Frame::applyXDoG(float sigma, float k)
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
            //Threshold DoG
            float dog_val = threshold(total1, total2, 0);

            //XDoG
            float p = 4.7;

            float g1, g2;
            g1 = total1 / 255;
            g2 = total2 / 255;

            float u = (1 + p) * (g1) - p * (g2); //Sharpened Image
            float t = 1;     //Threshold
            float e = 0.3;   //Level above which is white
            float phi = 0.01;//0.027; //Sharpness of transitions from black to white.

            if(u < e)
                t = tanh(phi * (u - e));

            dog_val = (1 - t) * g1 + t * (g1-g2);
            if(dog_val < 0.03)
                dog_val = 1;
            else
                dog_val = 0;
            dog_val*=255;
            
            dog_val = (unsigned char)fmax(fmin(255, dog_val), 0);

            new_img[pos + 0] =
                new_img[pos + 1] =
                new_img[pos + 2] = dog_val;
        }

    }
    delete data;
    delete k1;
    delete k2;
    data = new_img;
    
    return true;
}

/* applyColorXDoG(float sigma, float k)
 * 
 * PURPOSE:
 * Applies a XDoG filter to the image
 * XDoG is a DoG filter augmented with values for customization
 * and smoother color transitions.
 * Thresholds the result to either black or the original pixel color
 * 
 * INPUT:
 * float sigma: The size of Kernel A
 * float k: The scalar with which to augment the size of Kernel B
 * 
 * EXPECTATIONS:
 * It is preferable to have the image pre-filtered with a bilateral filter
 * k should be 1 or above, (1 will have no effect)
 *
 * NOTE: 
 * Larger sigma or k values will have a negative impact on performance
 * 
 * OUTPUT:
 * Applies XDoG edge finding to global image: char* data 
 * 
 */

bool Frame::applyColorXDoG(float sigma, float k)
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
            //Threshold DoG
            float dog_val = threshold(total1, total2, 0);

            //XDoG
            float p = 25;//4.7;

            float g1, g2;
            g1 = total1 / 255;
            g2 = total2 / 255;

            float u = (1 + p) * (g1) - p * (g2); //Sharpened Image
            float t = 1;                         //Threshold
            float e = 0.2;                       //Level above which is white
            float phi = 0.1;//0.01;//0.027;            //Sharpness of transitions from black to white.

            if(u < e)
                t = tanh(phi * (u - e));

            dog_val = (1 - t) * g1 + t * (g1 - g2);
            if(dog_val < 0.02){
                new_img[pos + 0] = data[pos + 0];
                new_img[pos + 1] = data[pos + 1];
                new_img[pos + 2] = data[pos + 2];
            } else {
                dog_val = 0;

                new_img[pos + 0] =
                    new_img[pos + 1] =
                    new_img[pos + 2] = dog_val;
            }
        }

    }
    delete data;
    delete k1;
    delete k2;
    data = new_img;
    
    return true;
}

/* applyColorDoG(float sigma, float k)
 * 
 * PURPOSE:
 * Applies a DoG filter to the image
  * Thresholds the result to either black or the original pixel color
 * 
 * INPUT:
 * float sigma: The size of Kernel A
 * float k: The scalar with which to augment the size of Kernel B
 * 
 * EXPECTATIONS:
 * It is preferable to have the image pre-filtered with a bilateral filter
 * k should be 1 or above, (1 will have no effect)
 *
 * NOTE: 
 * Larger sigma or k values will have a negative impact on performance
 * 
 * OUTPUT:
 * Applies DoG edge finding to global image: char* data 
 * 
 */

bool Frame::applyColorDoG(float sigma, float k)
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

            float dog_val = threshold(total1, total2, 0);
            float p = 15.7;
            float u = dog_val = (1 + p) * total1 - p * total2;
            float t = 1;
            float e = 200.5;
            float phi = 0.027;

            if(u < e)
                t =  tanh(phi * (u - e));

            if(t != 1){
                dog_val = (1 + t) * total1 - t * (total1 - total2);
            
                dog_val = (unsigned char)fmax(fmin(255, dog_val), 0);

                new_img[pos + 0] =
                    new_img[pos + 1] =
                    new_img[pos + 2] = dog_val;
            } else {
                new_img[pos + 0] = data[pos + 0];
                new_img[pos + 1] = data[pos + 1];
                new_img[pos + 2] = data[pos + 2];
            }
            

        }

    }
    delete data;
    delete k1;
    delete k2;
    data = new_img;
    
    return true;
}

float Frame::threshold(float a, float b, float t){
    float res;
    if((a - b) > t){
        res = 255;// - (a - b);
    } else {
        res = 0;
    }
    return res;
}

unsigned char Frame::get(int x, int y, int c)
{
    return (x < width && y < height && x >= 0 && c < channels && y >= 0 && c >= 0) ? data[y*width*channels + x*channels + c] : 0;
}

/* applyKuwahara(int a)
 * 
 * PURPOSE:
 * Applies a Kuwahara filter to the image
 * 
 * INPUT:
 * float a: width/height of each region explored
 * 
 * EXPECTATIONS:
 * n/a
 *
 * NOTE: 
 * n/a
 * 
 * OUTPUT:
 * Applies kuwahara filter to global image: char* data 
 * 
 */

bool Frame::applyKuwahara(int a)
{
    unsigned char *new_img = (unsigned char *)calloc(channels*width*height, sizeof(char));

    if(!new_img) {
        return false;
    }

	for (int iy=0 ; iy<height ; iy++) {
		for (int ix=0 ; ix<width ; ix++) {
            int pos = iy*width*channels + ix*channels;
            double stdd[4] = {0};
            double mins = 999999999999;
            int min = -1;

            double meanL[4] = {0};
            double meanR[4] = {0};
            double meanG[4] = {0};
            double meanB[4] = {0};

            // Over our four quadrants
            for(int q=0; q<4; q++) {
                int ky, h, kx, w;

                // One-pass standard deviation calculation
                double M = 0;
                double Q = 0;
                int count = 0;

                // Quadrant boundaries
                if(q < 2) {
                    ky = fmax(0, iy);
                    h = fmin(height, iy+a+1);
                }
                else {
                    ky = fmax(0, iy-a);
                    h = fmin(height, iy+1);
                }

                for(; ky<h; ky++) {
                    // Quadrant boundaries
                    if(q == 0 || q == 3) {
                        kx = fmax(0, ix);
                        w = fmin(width, ix+a+1);
                    }
                    else {
                        kx = fmax(0, ix-a);
                        w = fmin(width, ix+1);
                    }

                    for(; kx<w; kx++) {
                        // Luminance Value
                        int kpos = ky*width*channels + kx*channels;
                        double val = data[kpos + 0]*0.21 + data[kpos + 1]*0.71 + data[kpos + 2]*0.07;
                        Q = (count == 0) ? 0 : Q + (count-1)*(val-M)*(val-M)/count;
                        M = (count == 0) ? val : M + (val-M)/count;
                        meanL[q] += val;
                        meanR[q] += data[kpos + 0];
                        meanG[q] += data[kpos + 1];
                        meanB[q] += data[kpos + 2];
                        count++;
                    }
                }

                stdd[q] = Q/count;
                meanL[q] /= count;
                meanR[q] /= count;
                meanG[q] /= count;
                meanB[q] /= count;

                if(stdd[q] < mins) {
                    mins = stdd[q];
                    min = q;
                }
            }

            new_img[pos + 0] = (unsigned char)fmax(fmin(255, meanR[min]), 0);
            new_img[pos + 1] = (unsigned char)fmax(fmin(255, meanG[min]), 0);
            new_img[pos + 2] = (unsigned char)fmax(fmin(255, meanB[min]), 0);
		}
	}
    delete data;
    data = new_img;
    
    return true;
}

/* applyAnisotropicKuwahara(int radius, float q)
 * 
 * PURPOSE:
 * Applies an anistropic kuwahara filter to image
 * 
 * INPUT:
 * int radius: 
 * float q:
 * 
 * EXPECTATIONS:
 * n/a
 *
 * NOTE: 
 * DOES NOT CURRENTLY WORK!!
 * Much of the current code was heavily based on GLSL implementations in an attempt to have a working
 * filter by the due date. See Report for more information.
 * 
 * OUTPUT:
 * Applies Kuwahara filter to global image: char* data 
 * 
 */

bool Frame::applyAnisotropicKuwahara(int radius, float q)
{
    unsigned char *new_img = (unsigned char *)calloc((channels+1)*width*height, sizeof(char));
    unsigned char *tfm = (unsigned char *)calloc((channels+1)*width*height, sizeof(char));
    Frame K0123; K0123.open("../krnlx4_32n8.png");

    unsigned char pix[3];

    if(!new_img || !tfm) {
        return false;
    }

	for (int iy=0 ; iy<height ; iy++) {
		for (int ix=0 ; ix<width ; ix++) {
            int pos = iy*width*channels + ix*channels;
            int pos2 = iy*width*(channels+1) + ix*(channels+1);

            pix[0] = data[pos + 0];
            pix[1] = data[pos + 1];
            pix[2] = data[pos + 1];

            float lambda1 = 0.5 * (pix[1] + pix[0] + sqrt(pix[1]*pix[1] - 2.0*pix[0]*pix[1] + pix[0]*pix[0] + 4.0*pix[2]*pix[2]));
            float lambda2 = 0.5 * (pix[1] + pix[0] - sqrt(pix[1]*pix[1] - 2.0*pix[0]*pix[1] + pix[0]*pix[0] + 4.0*pix[2]*pix[2]));

            float v[2] = {lambda1 - pix[0], - pix[2]};

            if (LENGTH_2(v) > 0.0) {
                v[0] /= (float)LENGTH_2(v);
                v[1] /= (float)LENGTH_2(v);
            } else {
                v[0] = 0.0;
                v[1] = 1.0;
            }

            float phi = atan(v[1]/v[0]);
            float A = (lambda1 + lambda2 > 0.0) ? (lambda1 - lambda2) / (lambda1 + lambda2) : 0.0;

            tfm[pos2 + 0] = (unsigned char)fmax(fmin(255, v[0]), 0);
            tfm[pos2 + 1] = (unsigned char)fmax(fmin(255, v[1]), 0);
            tfm[pos2 + 2] = (unsigned char)fmax(fmin(255, phi), 0);
            tfm[pos2 + 3] = (unsigned char)fmax(fmin(255, A), 0);
        }
    }

	for (int iy=0 ; iy<height ; iy++) {
		for (int ix=0 ; ix<width ; ix++) {
            int pos = iy*width*channels + ix*channels;
            int pos2 = iy*width*(channels+1) + ix*(channels+1);

            //src_size = vec2(textureSize2D(src, 0));
            //vec2 uv = gl_FragCoord.xy / src_size;

            float m[8][4];
            float s[8][3];

            for(int ii=0; ii<8; ii++) {
                int jj = 0;
                for(jj=0; jj<4; jj++) {
                    m[ii][jj] = 0;
                    s[ii][jj] = 0;
                }
                m[ii][jj] = 0;
            }

            float t[4] = {tfm[pos2 + 0], tfm[pos2 + 1], tfm[pos2 + 2], tfm[pos2 + 3]};
            float c[4] = {data[pos + 0], data[pos + 1], data[pos + 2], data[pos + 3]};

            float a = radius;// * clamp((alpha + t.w) / alpha, 0.1, 2.0);
            float b = radius;// * clamp(alpha / (alpha + t.w), 0.1, 2.0);

            float cos_phi = cos(t[2]);
            float sin_phi = sin(t[2]);

            float R[4] = {cos_phi, -sin_phi, sin_phi, cos_phi};
            float S[4] = {0.5/a, 0.0, 0.0, 0.5/b};
            float SR[4] = MAT_MULT_22(S, R);

            int max_x = int(sqrt(a*a * cos_phi*cos_phi + b*b * sin_phi*sin_phi));
            int max_y = int(sqrt(a*a * sin_phi*sin_phi + b*b * cos_phi*cos_phi));

            //vec3 c = texture2D(src, uv).rgb;
            float w = K0123.get(K0123.getHeight()/2.0f, K0123.getWidth()/2.0f, 0);

            for (int k = 0; k < 8; ++k) {
                m[k][0] += c[0] * w;
                m[k][1] += c[1] * w;
                m[k][2] += c[2] * w;
                m[k][3] += w;
                s[k][0] += c[0] * c[0] * w;
                s[k][1] += c[1] * c[1] * w;
                s[k][2] += c[2] * c[2] * w;
            }

            for (int j = 0; j <= max_y; ++j) {
                for (int i = -max_x; i <= max_x; ++i) {
                    if ((j !=0) || (i > 0)) {
                        //float v[2] = {SR*i, SR*j};
                        float temp[2] = {i,j};
                        float v[2] = MAT_MULT_21(SR, temp);

                        if (DOT_2(v,v) <= 0.25) {
                            int kpos = pos + (i*width*channels + j*channels);
                            float c0[3] = {data[kpos + 0], data[kpos + 1], data[kpos + 2]};
                            kpos = pos - (i*width*channels + j*channels);
                            float c1[3] = {data[kpos + 0], data[kpos + 1], data[kpos + 2]};
                            //vec3 c0 = texture2D(src,uv + vec2(i,j)/src_size).rgb;
                            //vec3 c1 = texture2D(src,uv - vec2(i,j)/src_size).rgb;
                            //vec3 cc0 = c0 * c0;
                            //vec3 cc1 = c1 * c1;
                            float cc0[3] = {c0[0]*c0[0], c0[1]*c0[1], c0[2]*c0[2]};
                            float cc1[3] = {c1[0]*c1[0], c1[1]*c1[1], c1[2]*c1[2]};

                            int kix = (0.5+v[0])/K0123.getWidth();
                            int kiy = (0.5+v[1])/K0123.getWidth();
                            float w0123[4] = {K0123.get(kix, kiy, 0),K0123.get(kix, kiy, 1),K0123.get(kix, kiy, 2),K0123.get(kix, kiy, 3)};
                            //vec4 w0123 = texture2D(K0123, vec2(0.5, 0.5) + v);
                            for (int k = 0; k < 4; ++k) {
                                m[k][0] += c0[0]*w0123[k];
                                m[k][1] += c0[1]*w0123[k];
                                m[k][2] += c0[2]*w0123[k];
                                m[k][3] += w0123[k];
                                s[k][0] += cc0[0]*w0123[k];
                                s[k][1] += cc0[1]*w0123[k];
                                s[k][2] += cc0[2]*w0123[k];
                                //m[k] = {m[k][0] + c0[0]*w0123[k], m[k][1] + c0[1]*w0123[k], mk[2] + c0[2]*w0123[k], m[3] + w0123[k]};
                                //s[k] = {s[k][0] + cc0[0]*w0123[k], s[k][1] + cc0[1]*w0123[k], sk[2] + cc0[2]*w0123[k]};

                                //vec4(c0 * w0123[k], w0123[k]);
                                //s[k] += cc0 * w0123[k];
                            }
                            for (int k = 4; k < 8; ++k) {
                                m[k][0] += c1[0]*w0123[k-4];
                                m[k][1] += c1[1]*w0123[k-4];
                                m[k][2] += c1[2]*w0123[k-4];
                                m[k][3] += w0123[k-4];
                                s[k][0] += cc1[0]*w0123[k-4];
                                s[k][1] += cc1[1]*w0123[k-4];
                                s[k][2] += cc1[2]*w0123[k-4];
                                //m[k] = {m[k][0] + c1[0]*w0123[k-4], m[k][1] + c1[1]*w0123[k-4], mk[2] + c1[2]*w0123[k-4], m[3] + w0123[k-4]};
                                //s[k] = {s[k][0] + cc1[0]*w0123[k-4], s[k][1] + cc1[1]*w0123[k-4], sk[2] + cc1[2]*w0123[k-4]};

                                //m[k] += vec4(c1 * w0123[k-4], w0123[k-4]);
                                //s[k] += cc1 * w0123[k-4];
                            }

                            kix = (0.5-v[0])/K0123.getWidth();
                            kiy = (0.5-v[1])/K0123.getWidth();
                            float w4567[4] = {K0123.get(kix, kiy, 0),K0123.get(kix, kiy, 1),K0123.get(kix, kiy, 2),K0123.get(kix, kiy, 3)};
                            //vec4 w4567 = texture2D(K0123, vec2(0.5, 0.5) - v);
                            for (int k = 4; k < 8; ++k) {
                                m[k][0] += c0[0]*w4567[k-4];
                                m[k][1] += c0[1]*w4567[k-4];
                                m[k][2] += c0[2]*w4567[k-4];
                                m[k][3] += w4567[k-4];
                                s[k][0] += cc0[0]*w4567[k-4];
                                s[k][1] += cc0[1]*w4567[k-4];
                                s[k][2] += cc0[2]*w4567[k-4];
                                //m[k] = {m[k][0] + c0[0]*w4567[k-4], m[k][1] + c0[1]*w4567[k-4], mk[2] + c0[2]*w4567[k-4], m[3] + w4567[k-4]};
                                //s[k] = {s[k][0] + cc0[0]*w4567[k-4], s[k][1] + cc0[1]*w4567[k-4], sk[2] + cc0[2]*w4567[k-4]};


                                //m[k] += vec4(c0 * w4567[k-4], w4567[k-4]);
                                //s[k] += cc0 * w4567[k-4];
                            }
                            for (int k = 0; k < 4; ++k) {
                                m[k][0] += c1[0]*w4567[k];
                                m[k][1] += c1[1]*w4567[k];
                                m[k][2] += c1[2]*w4567[k];
                                m[k][3] += w4567[k];
                                s[k][0] += cc1[0]*w4567[k];
                                s[k][1] += cc1[1]*w4567[k];
                                s[k][2] += cc1[2]*w4567[k];

                                //m[k] += vec4(c1 * w4567[k], w4567[k]);
                                //s[k] += cc1 * w4567[k];
                            }
                        }
                    }
                }
            }

            //vec4 o = vec4(0.0);
            float o[4] = {0};
            for (int k = 0; k < 8; ++k) {
                m[k][0] /= m[k][3];
                m[k][1] /= m[k][3];
                m[k][2] /= m[k][3];
                s[k][0] = abs(s[k][0]/m[k][3] - m[k][0]*m[k][0]);
                s[k][1] = abs(s[k][1]/m[k][3] - m[k][1]*m[k][1]);
                s[k][2] = abs(s[k][2]/m[k][3] - m[k][2]*m[k][2]);
                
                //m[k].rgb /= m[k].w;
                //s[k] = abs(s[k] / m[k].w - m[k].rgb * m[k].rgb);

                //float sigma2 = sqrt(s[k].r) + sqrt(s[k].g) + sqrt(s[k].b);
                //float w = 1.0 / (1.0 + pow(255.0 * sigma2, q));
                float sigma2 = sqrt(s[k][0]) + sqrt(s[k][1]) + sqrt(s[k][2]);
                float w = 1.0f / (1.0f + pow(255.0f * sigma2, q));

                //o += vec4(m[k].rgb * w, w);
                o[0] += m[k][0]*w;
                o[1] += m[k][1]*w;
                o[2] += m[k][2]*w;
                o[3] += w;
            }
            //printf("o = %f %f %f %f\n", o[0], o[1], o[2], o[3]);

            //gl_FragColor = vec4(o.rgb / o.w, 1.0);
            new_img[pos + 0] = (unsigned char)fmax(fmin(255, o[0]/o[3]), 0);
            new_img[pos + 1] = (unsigned char)fmax(fmin(255, o[1]/o[3]), 0);
            new_img[pos + 2] = (unsigned char)fmax(fmin(255, o[2]/o[3]), 0);
        }
    }

    delete tfm;
    delete data;
    data = new_img;
    
    return true;
}

