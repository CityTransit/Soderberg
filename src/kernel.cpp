#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "kernel.h"


Kernel* Kernel::generateGaussian(int sigma, int std_dev, int distance){
    float total = 0;
    float *kernel = (float *)malloc(sizeof(float) * distance * distance);

    for(int i = 0; i < distance; i++){
        for(int j = 0; j < distance; j++){
            kernel[distance * j + i] = sigma;
        }
    }

    for(int i = 0; i < distance; i++){
        for(int j = 0; j < distance; j++){
            float g =  exp(-1 / (2* pow(std_dev, 2)) * (pow(i, 2) + pow(j, 2) ));
            total += g;
            kernel[distance * j + i] = g;
        }
    }
    Kernel *k = new Kernel(distance, distance, kernel, total); 
    return k;
}

Kernel* Kernel::generateBilateral(float s, float r, int size){
    float total = 0;
    float *kernel = (float *)malloc(sizeof(float) * size * size);

    for(int x = 0; x < size; x++){
        for(int y = 0; y < size; y++){
            float sum = 0;
            
            for(int i = 0; i < size; i++){
                for(int j = 0; j < size; j++){
                    float delta = sqrt(pow((x - i), 2) + pow((y - j), 2));
                    float abs = fabs(kernel[size * x + y] - kernel[size * i + j]);
                    sum+=gaussian(delta, s) * gaussian(abs, r);
                }
            }

            sum/=size*size;
            total += sum;
            kernel[size * x + y] = sum;
            printf("%f ", sum);
        }
        printf("\n");
    }

    Kernel *k = new Kernel(size, size, kernel, total); 
    return k;
}

float Kernel::gaussian(float n, float sigma){
    return exp(-1 / (2* pow(sigma, 2)) * ( pow(n, 2) ));
}
