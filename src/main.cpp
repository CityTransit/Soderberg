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

    f.open("test_image.png");

    k = Kernel::generateGaussian(4, 5);
    f.applyKernel(k);

    f.save("test_out.png", "butts");

    return 0;
}
