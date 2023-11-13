#include "stdio.h"
#include "stdlib.h"
#include "wsart_support.h"
// #include "png.h"

/*
    Program flow                                ----------
                                                V        |
    Read image --> make sinogram --> wsart --> wbp --> wfp --> wbp --> out

        A
        |
        |
   On this step

*/

int main(int argc, char const *argv[]) {
    
    const char *in_file_name;
    //const char *sin_file_name = "sinogram.png";
    const char *out_file_name;
    
    cudaError_t cuda_ret = cudaErrorUnknown;
    //printf("got here");
   // exit(0);
    
    if(argc == 1) {
        printf("\nERROR: Missing Input file and output file arguments\n\n");
        printf("\tExpected 2 arguments, received 0\n\n");
        return 1;
    } else if(argc == 2) {
        printf("\nERROR: Missing Input file or output file argument\n\n");
        printf("\tExpected 2 arguments, received 1\n\n");
        return 2;
    } else if(argc == 3) {
        in_file_name = argv[1];
        out_file_name = argv[2];
    } else {
        printf("\n    Invalid input parameters!"
               "\n    Usage: ./wsart <input_file> <output_file>"
               "\n\n");
        exit(0);
    }

    img_t inimg;
    img_t outimg;
    //printf("got here");
    //exit(0);
    read_tiff(in_file_name, &inimg);

    copyimg(&inimg, &outimg);

    write_tiff(out_file_name, &outimg);

    return 0;
}
