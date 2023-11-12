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
    

    FILE *in_file;
    FILE *sin_file;
    FILE *out_file;
    const char *in_file_name;
    const char *sin_file_name = "sinogram.png";
    const char *out_file_name;
    png_t *in_png;
    png_t *sin_png;
    png_t *out_png;

    cudaError_t cuda_ret;
    if(argc == 1) {
        printf("\nERROR: Missing Input file and output file arguments\n");
        printf("Expected 2 arguments, received 0");
        return 1;
    } else if(argc == 2) {
        printf("\nERROR: Missing Input file or output file argument\n");
        printf("Expected 2 arguments, received 1");
        return 2;
    } else if(argc == 3) {
        in_file_name = argv[1];
        out_file_name = argv[2];
    } else {
        printf("\n    Invalid input parameters!"
               "\n    Usage: ./wsart <input_file> <output_file>"
               "\n");
        exit(0);
    }

    in_file = fopen(in_file_name, "rb");
    //fopen sinfile
    //fopen outfile


    fclose(in_file);
    //fclose sinfile
    //fclose outfile

    return 0;
}
