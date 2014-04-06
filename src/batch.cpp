#ifdef USE_OPEN_MP
#include <omp.h>
#endif

#include "png.h"
#include "frame.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "kernel.h"

int main(int argc, char *argv[])
{
    if(argc < 6) {
        printf("ERROR: Invalid number of arguments.\n\tUsage: %s start end prefix in_folder out_folder\n", argv[0]);
        return 0;
    }

    char outname[100];
    const int start = atoi(argv[1]);
    const int end = atoi(argv[2]);
    const char *prefix = argv[3];
    const char *in_folder = argv[4];
    const char *out_folder = argv[5];

    printf ("Opening files...\n");

#ifdef USE_OPEN_MP
    #pragma omp parallel for private(outname)
#endif
    for(int i=start; i<end; i++) {

            Frame f;

            snprintf(outname, 100, "%s/%s%05d.png", in_folder, prefix, i);
            f.open(outname);

            printf("Processing file %s...\n", outname);
            //f.applyKuwahara(3);
            f.applyColorDoG(4, 1.6);

            snprintf(outname, 100, "%s/%s%05d.png", out_folder, prefix, i);
            f.save(outname, outname);
    }

    return 0;
}

