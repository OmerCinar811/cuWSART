#include "png.h"

png_infop info_ptr;
png_bytepp row_pointers;

void read_png(char*);
void write_png(char*);