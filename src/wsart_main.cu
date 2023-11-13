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
    //FILE *sin_file;
    FILE *out_file;
    const char *in_file_name;
    //const char *sin_file_name = "sinogram.png";
    const char *out_file_name;
    png_t in_png;
    //png_t sin_png;
    png_t out_png;

    //cudaError_t cuda_ret;
    
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

    in_file = fopen(in_file_name, "rb");
    out_file = fopen(out_file_name, "wb");
    //fopen sinfile

    // read_png(in_file, &in_png);
    // memcpy(&out_png, &in_png, sizeof(in_png));
    // write_png(out_file, &out_png);
 

    fclose(in_file);
    fclose(out_file);
    //fclose sinfile

    return 0;
}
