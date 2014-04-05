/*
 *
 * Code ruthlessly stolen from:
 *  - http://www.labbookpages.co.uk/software/imgProc/libPNG.html
 *  - http://www.piko3d.net/tutorials/libpng-tutorial-loading-png-files-from-streams/#Loading
*/

#include <stdint.h>
#include <math.h>
#include <omp.h>

#include <dirent.h>
#include "png.h"
#include "frame.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "kernel.h"

int main(int argc, char *argv)
{
    if(argc < 4) {
        printf("ERROR: Invalid number of arguments.\n\tUsage: %s num_frames in_folder out_folder\n", argv[0]);
    }

    struct dirent *ent;
    DIR *dir;
    char outname[50];
    const int num_frames = atoi(argv[1]);
    const char *in_folder = argv[2];
    const char *out_folder = argv[3];

    printf ("Opening files...");
    if ((dir = opendir (in_folder)) != NULL) {
        readdir(dir); // .
        readdir(dir); // ..

        #pragma omp parallel for private(ent)
        for(int i=0; i<num_frames; i++) {
            #pragma omp critical
            {
                 ent = readdir(dir);
            }
    
            if ( ent != NULL) {
                Frame f;
                snprintf(outname, 100, "%s/%s", in_folder, ent->d_name);
                f.open(outname);
                printf("Processing file %s...\n", outname);

                f.applyKuwahara(3);

                snprintf(outname, 100, "%s/%s", out_folder, ent->d_name);
                f.save(outname, outname);
            }
        }
        closedir (dir);
    } else {
        fprintf(stderr, "Error: Unable to open input directory.\n");
        return EXIT_FAILURE;
    }

    return 0;
}

