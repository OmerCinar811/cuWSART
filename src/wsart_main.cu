#include "stdio.h"
#include "stdlib.h"
#include "wsart_support.h"
// #include "png.h"

/*
    Program flow                                ----------
                                                V        |
    Read image --> make sinogram --> wsart --> wbp --> wfp --> wbp --> out

*/

/*
int main() {



}
*/

int main(int argc, char *argv[])
{
    read_png(argv[1]);
    write_png(argv[2]);
    return 0;
}