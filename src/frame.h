#pragma once

#include "png.h"
#include "kernel.h"
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>

#define PNGSIGSIZE 8

class Frame
{
public:
    Frame() {}
    Frame(Frame *copy);
    ~Frame();

    bool open(const char *fname);
    bool save(const char *filename, char *title);                           // Save image to file

    bool flip();

    bool applyKernel(Kernel *k);
    bool applyTwoToneKernel(Kernel *k);
    bool applyDoG(float s, float k);
    bool applyXDoG(float s, float k);
    bool applyColorXDoG(float s, float k);
    bool applyColorDoG(float s, float k);

    bool applyBilateral(float s, float r);
    bool applyTemporalBilateral(float s, float r, float t, Frame *f, int fn);

    bool applyKuwahara(int a);
    bool applyAnisotropicKuwahara(int radius, float q);

    float threshold(float a, float b, float t);
    static float gaussian(float n, float sigma);

    png_uint_32 getWidth() {return width;}
    png_uint_32 getHeight() {return height;}
private:
    static void userReadData(png_structp pngPtr, png_bytep data, png_size_t length);
    static bool validate(std::istream &src);
    unsigned char get(int x, int y, int c);

    png_uint_32 width;
    png_uint_32 height;
    png_uint_32 bitdepth;
    png_uint_32 channels;
    png_uint_32 colourtype;

    png_structp pngPtr;
    png_bytep *rowPtrs;
    png_infop infoPtr;
    unsigned char *data;
}; 
