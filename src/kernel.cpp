#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.h"

Kernel *Kernel::generateGaussian(int sigma, int std_dev, int distance)
{
    float total = 0;
    float *kernel = (float *)malloc(sizeof(float) * distance * distance);

    for(int i = 0; i < distance; i++) {
        for(int j = 0; j < distance; j++) {
            float g =  exp(-1 / (2* pow(std_dev, 2)) * (pow(i, 2) + pow(j, 2) ));
            total += g;
            kernel[distance * j + i] = g;
        }
    }

    return new Kernel(distance, distance, kernel, total); 
}

Kernel *Kernel::generateSharpen() 
{
    int width = 3;
    int height = 3;
    float *values = (float *)calloc(width*height, sizeof(float));

    values[1] = -1;
    values[3] = -1;
    values[5] = -1;
    values[7] = -1;
    values[4] = 5;

    return new Kernel(width, height, values);
}

Kernel *Kernel::generateEdge()
{
    int width = 3;
    int height = 3;
    float *values = (float *)malloc(width*height*sizeof(float));


    for(int i=0; i<9; i++) values[i] = -1;
    values[4] = 8;

    return new Kernel(width, height, values);
}

Kernel *Kernel::generateEmboss()
{
    int width = 3;
    int height = 3;
    float *values = (float *)malloc(width*height*sizeof(float));


    values[0] = -2;
    values[1] = -1;
    values[2] = 0;
    values[3] = -1;
    values[4] = 1;
    values[5] = 1;
    values[6] = 0;
    values[7] = 1;
    values[8] = 2;

    return new Kernel(width, height, values);
}

