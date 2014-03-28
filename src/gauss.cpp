#include <math.h>

#include "gauss.h"


Gauss::Gauss(float r){
    for(int i = 0; i < 256; i++){
        g[i] = gaussian(i, r);
    }
}

float Gauss::getGauss(int x){
    return g[x];
}

float Gauss::gaussian(float x, float sigma){
    return exp(- x*x / (2*sigma*sigma) ) / (2*M_PI*sigma);
}
