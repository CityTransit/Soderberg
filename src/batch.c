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

int main()
{
    Frame f;
    //Kernel *k = Kernel::generateGaussian(100, 5, 10);
    Kernel *k = Kernel::generateLoG(5, 10);
    DIR *dir;
    struct dirent *ent;
    char outname[50];

    /* print all the files and directories within directory */
    if ((dir = opendir ("../input/")) != NULL) {
        readdir(dir);
        readdir(dir);
        while ((ent = readdir(dir)) != NULL) {
              printf ("Opening %s...\n", ent->d_name);
              snprintf(outname, 50, "../input/%s", ent->d_name);
              f.open(outname);
              f.applyKernel(k);
              snprintf(outname, 50, "../output/%s", ent->d_name);
              f.save(outname, outname);
        }
        closedir (dir);
    } else {
    /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
    }

    return 0;
}

