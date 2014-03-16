/*
 *
 * Code ruthlessly stolen from:
 *  - http://www.labbookpages.co.uk/software/imgProc/libPNG.html
 *  - http://www.piko3d.net/tutorials/libpng-tutorial-loading-png-files-from-streams/#Loading
*/


#include "png.h"
#include "frame.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "kernel.h"

int main()
{
    Frame f;
    Kernel *k;

    float *d = (float *)malloc(9*sizeof(float));

    for(int i=0; i<9; i++) d[i] = 1;

    //k = new Kernel(3, 3, d, 9);

    //k = Kernel::generateGaussian(100, 5, 10);
    k = Kernel::generateBilateral(4, 0.05, 10);
    
    f.open("test_image.png");

    //f.flip();
    //f.applyKernel(k);
    f.applyBilateral(4, 0.2, 10);


    f.save("test_out.png", "butts");

    return 0;
}
