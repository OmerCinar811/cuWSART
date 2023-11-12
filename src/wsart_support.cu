#include "wsart_support.h"

// png_infop info_ptr;
// png_bytepp row_pointers;

// void read_png(char *file_name) {
//     FILE *fp = fopen(file_name, "rb");
//     png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//     info_ptr = png_create_info_struct(png_ptr);
//     png_init_io(png_ptr, fp);
//     png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
//     row_pointers = png_get_rows(png_ptr, info_ptr);
//     png_destroy_read_struct(&png_ptr, NULL, NULL);
//     fclose(fp);
// }

// void write_png(char *file_name) {
//     FILE *fp = fopen(file_name, "wb");
//     png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//     png_init_io(png_ptr, fp);
//     png_set_rows(png_ptr, info_ptr, row_pointers);
//     png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
//     png_destroy_write_struct(&png_ptr, &info_ptr);
//     fclose(fp);
// }

/**
 * @brief Reads in a png
 *        Assumes that file is a png
 * 
 * @param file_name string for file name
 */
void read_png(char *file_name) {

    FILE *fp = fopen(file_name, "rb"); // open file

    if(!fp) { //check if file was able to be opened
        printf("Could not open file\n");
        exit(FILE_ERR);
    } 

    // Creating png read struct and checking if it was created
    png_structp p_png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if(!p_png) {
        printf("Could not create png read struct pointer\n");
        exit(PNG_STRUCT_ERR);
    }

    // Creating png info structure and checking if it was created
    png_infop p_info = png_create_info_struct(p_png);

    if(!p_info) {
        png_destroy_read_struct(&p_png, (png_infopp)NULL, (png_infopp)NULL); // free any mem of png structs
        printf("Could not create png info struct pointer\n");
        exit(PNG_INF_STR_ERR);
    } 

    // set file pointer to png struct
    png_init_io(p_png, fp);
    png_set_sig_bytes(p_png, 8); //tells library to ignore first 8 bytes (aka magic number)


    // read info from the png, including the info before the data and the header
    png_read_info(p_png, p_info);
    uint32_t width, height;
    int bit_depth, color_type, interlace_type, compression_type, filter_method;
    
    // Get image header data
    png_get_IHDR(p_png, p_info,  &width, &height, &bit_depth, &color_type, &interlace_type, NULL, NULL);


    // Read the actual png image data
    png_bytepp row_pointers = (png_bytepp)png_malloc(p_png, sizeof(png_bytepp)*height);
    for(int i = 0; i < height; i++) {
        row_pointers[i] = (png_bytep)png_malloc(p_png, width * 1);
    }

     

} // read_png