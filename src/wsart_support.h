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
void write_tiff(const char*, img_t*);

#endif