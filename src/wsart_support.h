#ifndef _WSART_SUPPORT_H
#define _WSART_SUPPORT_H

#include "stdint.h"
#include "tiffio.h"

typedef struct {
    uint64_t xmax;
    uint64_t ymax;
    float **data;
} img_t;


void read_tiff(const char*, img_t*);
int write_tiff(const char*, img_t*);
void copyimg(img_t*, img_t*);


#endif