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

<<<<<<< HEAD
    f.open("test_image.png");

    k = Kernel::generateGaussian(4, 5);
    f.applyKernel(k);
=======
    //float *d = (float *)malloc(9*sizeof(float));

    //for(int i=0; i<9; i++) d[i] = 1;

    //k = new Kernel(3, 3, d, 9);

    //k = Kernel::generateGaussian(100, 5, 10);

    k = Kernel::generateBilateral(4, 0.05, 10);

    //k = Kernel::generateEdge();
    
    f.open("test_image.png");

    //f.flip();
    //f.applyKernel(k);

    f.applyBilateral(4, 0.2, 10);

    //k = Kernel::generateSharpen();
    //k = Kernel::generateEmboss();
    //k = Kernel::generateGaussian(100, 5, 10);
    //f.applyKernel(k);
>>>>>>> 0b1d7577b92defdb542ea8568f40d2c21aed8525

    f.save("test_out.png", "butts");

    return 0;
}
