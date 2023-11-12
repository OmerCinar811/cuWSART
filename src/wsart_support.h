#ifndef _WSART_SUPPORT_H
#define _WSART_SUPPORT_H

#include "png.h"
#include "stdint.h"


#define FILE_ERR        1
#define PNG_STRUCT_ERR  2
#define PNG_INF_STR_ERR 3

typedef struct {
    png_structp p_png;
    png_infop p_info;
} png_t;


// void read_png(char*);
// void write_png(char*);



#endif