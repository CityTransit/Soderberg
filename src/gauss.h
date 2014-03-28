#pragma once

class Gauss
{
 public:
    Gauss(float r);
    ~Gauss();
    float g[256];
    float gaussian(float x, float sigma);
    float getGauss(int x);
};
