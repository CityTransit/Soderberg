#pragma once

class Kernel
{
public:
    Kernel(int w, int h, float *k, float n = 1) { width = w; height=h; this->k = k; norm = n; }
    ~Kernel();
    int get_height() { return height; }
    int get_width() { return width; }
    int get_norm() { return norm; }
    float get(int i) { return k[i]; }
    static Kernel *generateGaussian(float sigma, int distance);

    static Kernel* generateGaussian(int sigma, int std, int distance);
    static Kernel* generateBilateral(float s, float r, int size);
    static float gaussian(float n, float sigma);
    static Kernel *generateSharpen();
    static Kernel *generateEdge();
    static Kernel *generateEmboss();
    
private:
    int width;
    int height;
    int norm;
    float *k;
};
