/*
 *
 * Code ruthlessly stolen from:
 *  - http://www.labbookpages.co.uk/software/imgProc/libPNG.html
 *  - http://www.piko3d.net/tutorials/libpng-tutorial-loading-png-files-from-streams/#Loading
*/


#include <dirent.h>
#include "png.h"
#include "frame.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "kernel.h"

int main(int argc, char* argv[])
{
    if(argc < 3) exit(1);

    Frame f;
    //Kernel *kg = Kernel::generateGaussian(10);
    //Kernel *kg2 = Kernel::generateGaussian(12 * 2);
    //Kernel *k = Kernel::generateLoG(0.8, 5);
    
    f.open(argv[1]);
    //f.applyTwoToneKernel(k);
    //f.applyKernel(kg);
    //f.applyBilateral(32, 0.2);
    //char* test_file_name = "test_bilateral.png";
    //f.save(test_file_name, test_file_name);

    f.applyColorDoG(4, 1.6);
    //    f.applyColorXDoG(8, 1.6);
    f.save(argv[2], argv[2]);


    return 0;
}

