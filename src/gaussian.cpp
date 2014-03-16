#include "gaussian.h"

Kernel* Gaussian::generateKernel(int std_dev, int distance){
    float total = 0;
    float *kernel = (float *)malloc(sizeof(float) * distance * distance);

    for(int i = 0; i < distance; i++){
        for(int j = 0; j < distance; j++){
            kernel[distance * j + i] = std_dev;
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
