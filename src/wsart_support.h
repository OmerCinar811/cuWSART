#ifndef _WSART_SUPPORT_H
#define _WSART_SUPPORT_H

#include "stdint.h"
#include "tiffio.h"
#include "math.h"
#include "math_constants.h"

#define SQR(X)		((X)*(X))
#define FLOATT      1
#define RGBTRPT     2
#define ARRWID      512
#define PHI		(phi*0.017453293)

typedef struct {
    uint64_t xmax;
    uint64_t ymax;
    float **data;
} img_t;

typedef struct {
    unsigned char red;
    unsigned char green;
    unsigned char blue;
} rgbtriplet;

/* Methods for reading and writing images */
int TIFFReadContigStripData(TIFF*, int, int, void*);
void read_tiff(const char*, img_t*);
int write_tiff(const char*, img_t*);
void copyimg(img_t*, img_t*);

/* Methods to generate sinograms */
void make_sinogram(img_t*, int, int, int);
void allocate_image(img_t*, int, int, int);
double ireadbuf(img_t*, double, double);
double readbuf_flt (img_t*, int, int);
void handback(img_t*, img_t*);
void freebuf(img_t*);


#endif