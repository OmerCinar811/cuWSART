#ifndef _WSART_SUPPORT_H
#define _WSART_SUPPORT_H

#include "png.h"
#include "stdint.h"


#define FILE_ERR                    0x11
#define PNG_STRUCT_ERR              0x12
#define PNG_INF_STR_ERR             0x13
#define PNG_ENDINF_STR_ERR          0x14
#define PNG_STRUCT_W_ERR            0x15

typedef struct {
    png_structp p_png;
    png_structp p_png_w;
    png_infop p_info;
    png_infop p_endinfo;
    uint8_t** row_pointers;
} png_t;

png_t create_png_wrapper(FILE*);
void destroy_png_wrapper(png_t*);
void read_png(png_t*);


// void read_png(char*);
// void write_png(char*);



#endif