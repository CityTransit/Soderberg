#include <math.h>
#include <stdlib.h>
#include "kernel.h"

class Gaussian
{
public:
    Kernel* generateKernel(int std_dev, int distance);
private:
    float* kernel;
};
