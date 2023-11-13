#include "wsart_support.h"

//png_infop info_ptr;
//png_bytepp row_pointers;

void read_png(FILE *fp, png_t* png) {
    //FILE *fp = fopen(file_name, "rb");
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png->info_ptr = png_create_info_struct(png_ptr);
    
    png_init_io(png_ptr, fp);
    png_read_png(png_ptr, png->info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    png->row_pointers = png_get_rows(png_ptr, png->info_ptr);
    png->height = png_get_image_height(png_ptr, png->info_ptr);
    png->width = png_get_image_width(png_ptr, png->info_ptr);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    //fclose(fp);
}

void write_png(FILE *fp, png_t* png) {
    //FILE *fp = fopen(file_name, "wb");
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_init_io(png_ptr, fp);
    png_set_rows(png_ptr, png->info_ptr, png->row_pointers);
    png_write_png(png_ptr, png->info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    png_destroy_write_struct(&png_ptr, &(png->info_ptr));
    //fclose(fp);
}