#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.h"

Kernel *Kernel::generateGaussian(int distance)
{
    float *kernel = (float *)malloc(sizeof(float) * distance * distance);
    float norm = 0;

    for(int i = 0; i < distance; i++) {
        float y = -(distance-1)/2 + i;
        for(int j = 0; j < distance; j++) {
            float x = -(distance-1)/2 + j;
            float g =  exp(-(x*x + y*y) / (2*distance*distance));
            kernel[distance * i + j] = g;
            norm += g;
        }
    }

    return new Kernel(distance, distance, kernel, norm); 
}

Kernel *Kernel::generateLoG(float sigma, int distance)
{
    float *kernel = (float *)malloc(sizeof(float) * distance * distance);

    int pos;
    for(int i = 0; i < distance; i++) {
        float y = -(distance-1)/2 + i;
        for(int j = 0; j < distance; j++) {
            float x = -(distance-1)/2 + j;
            pos = distance * i + j;
            
            float log = exp(-(x*x + y*y) / (2*sigma*sigma));
            log = (1 - (x*x + y*y) / (2*sigma*sigma)) * log;
            log = -( 1/(M_PI * pow(sigma, 4))) * log;

            kernel[pos] = log;
            printf("%f ", log);

        }
        printf("\n");
    }

    return new Kernel(distance, distance, kernel); 
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

