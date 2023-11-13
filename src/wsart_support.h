#ifndef _WSART_SUPPORT_H
#define _WSART_SUPPORT_H

#include "png.h"
#include "stdint.h"

// Error definitions
#define FILE_ERR                    0x11
#define PNG_STRUCT_ERR              0x12
#define PNG_INF_STR_ERR             0x13
#define PNG_ENDINF_STR_ERR          0x14
#define PNG_STRUCT_W_ERR            0x15

typedef struct {
    png_infop info_ptr;
    uint32_t width, height;
    uint8_t** row_pointers;
} png_t;


// void read_png(char*);
// void write_png(char*);

void read_png(FILE*, png_t*);
void write_png(FILE*, png_t*);

#endif