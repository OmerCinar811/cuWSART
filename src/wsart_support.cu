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
 * @brief Create a png wrapper object
 * 
 * @param fp 
 * @return png_t 
 */
png_t create_png_wrapper(FILE *fp) {
    png_t out;

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

    png_structp p_png_w = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!p_png_w) {
        printf("Could not create png write struct pointer\n");
        exit(PNG_STRUCT_W_ERR);
    }

    // Creating png info structure and checking if it was created
    png_infop p_info = png_create_info_struct(p_png);
    if(!p_info) {
        png_destroy_read_struct(&p_png, (png_infopp)NULL, (png_infopp)NULL); // free any mem of png structs
        printf("Could not create png info struct pointer\n");
        exit(PNG_INF_STR_ERR);
    } 

    png_infop p_endinfo = png_create_info_struct(p_png);
    if (!p_endinfo) {
        png_destroy_read_struct(&p_png, &p_info, (png_infopp)NULL);
        exit(PNG_ENDINF_STR_ERR);
    }
    out.p_png = p_png;
    out.p_info = p_info;
    out.p_endinfo = p_endinfo;
    // set file pointer to png struct
    png_init_io(p_png, fp);
    png_init_io(p_png_w, fp);

    return out;
} // create

void destroy_png_wrapper(png_t *png) {
    png_destroy_read_struct(&(png->p_png), NULL, (png_infopp)NULL);
    png_destroy_write_struct(&(png->p_png_w), NULL);
} // destroy

/**
 * @brief Reads in a png
 *        Assumes that file is a png
 * 
 * @param file_name string for file name
 */
void read_png(png_t *png) {

    png_set_sig_bytes(png->p_png, 8); //tells library to ignore first 8 bytes (aka magic number)
    // read info from the png, including the info before the data and the header
    png_read_info(png->p_png, png->p_info);
    uint32_t width = 0;
    uint32_t height = 0;
    int bit_depth = 0;
    int color_type = 0;
    int interlace_type = 0;
    int compression_type = 0;
    int filter_method = 0;
    
    // Get image header data
    png_get_IHDR(png->p_png, png->p_info,  &width, &height, &bit_depth, &color_type, &interlace_type, NULL, NULL);
    
    // assign various information about png to struct
    png->height = height;
    png->width = width;
    png->hdr.bit_depth = bit_depth;
    png->hdr.color_type = color_type;
    png->hdr.interlace_type = interlace_type;
    png->hdr.compression_type = compression_type;
    png->hdr.filter_method = filter_method;

    // Read the actual png image data
    uint8_t **row_pointers = (uint8_t**)png_malloc(png->p_png, sizeof(png_bytepp)*height);
    for(int i = 0; i < height; i++) {
        row_pointers[i] = (uint8_t*)png_malloc(png->p_png, width * 1);
    }

    //png_set_rows(png->p_png, png->p_info, row_pointers);
    png_read_image(png->p_png, row_pointers);
    png_read_end(png->p_png, png->p_endinfo);

} // read_png

void write_png(png_t *png) {

    png_set_IHDR(png->p_png_w,
                 png->p_info,
                 png->width,
                 png->height,
                 png->hdr.bit_depth,
                 png->hdr.color_type,
                 png->hdr.interlace_type,
                 png->hdr.compression_type,
                 png->hdr.filter_method);

    png_set_rows(png->p_png_w, png->p_info, png->row_pointers);
    png_write_png(png->p_png_w, png->p_info, PNG_TRANSFORM_IDENTITY, NULL);
    

} // write_png