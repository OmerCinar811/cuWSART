#include "wsart_support.h"

int TIFFReadContigStripData(TIFF* tif, int elemsize, int bpp, void* dataptr) {
    unsigned char *buf;
    tsize_t scanlinesize = TIFFScanlineSize(tif);
    long tiffstripsize;
    unsigned char* p;
    unsigned short* ps;
    char s[255];
    int i,j;
    uint32 row, h, w, ofs;
    uint32 rowsperstrip = (uint32)-1;

	p = (unsigned char*)dataptr; ps= (unsigned short*)dataptr;

	tiffstripsize = TIFFStripSize(tif); /* usually scan line (=row) size * rows per strip */
	buf = (unsigned char *)_TIFFmalloc(tiffstripsize);
	if (buf) 
	{

		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

		printf("Image length %d, image width (height) %d\n",(int)h,(int)w);
		printf("Scan line size: %d, strip size %li, rows per strip %d\n",
					(int)scanlinesize,tiffstripsize,rowsperstrip);
		if (bpp==12)
		{
			printf("Performing expansion from 12 to 16 bits/sample\n");
		}
		else if (bpp==1)	/* 1 bit per pixel (fax) - expand */
		{
			printf("Performing expansion from 1 to 8 bits/sample\n");
		}
		else if ((bpp==16)  && (elemsize==6))
		{
			printf("Performing conversion from 48-bit RGB to Z stack\n");
		}
		

		ofs=0;
		for (row = 0; row < h; row += rowsperstrip) 
		{
			/* nrow is the number of rows to be read. At the end of the image we have less of them */
			uint32 nrows = (row+rowsperstrip > h ? h-row : rowsperstrip);

			/* Read the next strip, decode compression if necessary */
			tstrip_t strip = TIFFComputeStrip(tif, row, 0);
			if (TIFFReadEncodedStrip(tif, strip, buf, nrows*scanlinesize) < 0) 
			{
				printf("Error reading strip data");
				return -1;
			} 
		
			if (bpp==12)		/* 12 bit per pixel, must expand */
			{
				ofs=0;
				for (i=0; i<nrows*scanlinesize; i+=3)	/* one element is 3 bytes, thus the increment */
				{
					if (ofs < h*w-1)
					{
					
						/* 3 bytes in buffer (xx-yy-zz) become 2 shorts (0xxy-0yzz) */
					
						// dp (3, "Src row %5d  Dest offset %5d   line %5d   %2X%2X%2X to %4X-%4X\n",
						// 		row, ofs, i,
						// 		buf[i],buf[i+1],buf[i+2],
						// 		((unsigned short)buf[i] << 4) + (buf[i+1] >> 4),
						// 		(((unsigned short)buf[i+1] << 8) & 0x0f00) + buf[i+2] );
					
						ps[w*row+ofs  ] = ((unsigned short)buf[i] << 4) + (buf[i+1] >> 4);
						ps[w*row+ofs+1] = (((unsigned short)buf[i+1] << 8) & 0x0f00) + buf[i+2] ;
						ofs += 2;
					}
				}
			}
			else if (bpp==1)	/* 1 bit per pixel (fax) - expand */
			{
				ofs=0;
				for (i=0; i<nrows*scanlinesize; i++)
				{
					if (ofs < h*w-1)
					{
					
						/* expand each bit of buf into one byte of p. */

						for (j=0; j<8; j++)
							p[w*row+ofs+j  ] = (buf[i] >> (7-j)) & 0x01;	/* Big-endian version? */
						ofs += 8;
					}
				}
			
			}
			else if ((bpp=16) && (elemsize==6))	/* 48-bit color RGB */
			{
				ofs=0;
				for (i=0; i<nrows*scanlinesize; i+=6)	/* one element is 3 shorts, thus the increment */
				{
					if (ofs < h*w-1)
					{
					
						/* 3 shorts in buffer are (rrr-ggg-bbb) */

						ps[w*row+ofs      ] = ((unsigned short)(buf[i+1]<<8) + buf[i+0]);
						ps[w*row+ofs+  h*w] = ((unsigned short)(buf[i+3]<<8) + buf[i+2]);
						ps[w*row+ofs+2*h*w] = ((unsigned short)(buf[i+5]<<8) + buf[i+4]);
						ofs += 1;
					}
				}
			
			}
			else			/* No expansion, copy verbatim */
			{
				memcpy (p+elemsize*row*w, buf, nrows*scanlinesize);
			}
			
		}
		_TIFFfree(buf);
	}
	else
		return -1;
		
	return 0;
}

void read_tiff(const char *filename, img_t *tiff) {

    TIFF *tiff_handle;
    uint32_t xsize, ysize;
    int bytes;
    uint16_t bitspersample, samplesperpixel, planarconfig, photometric, sampleformat;
    char s[256];
    double *pd; float *pf; long idx;

    tiff_handle = TIFFOpen(filename, "r");
    if(!tiff_handle) {
        printf("Unable to open tiff\n\n");
        exit(0x10);
    }

    TIFFGetField (tiff_handle, TIFFTAG_IMAGEWIDTH, &xsize);
	TIFFGetField (tiff_handle, TIFFTAG_IMAGELENGTH, &ysize);
	TIFFGetField (tiff_handle, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	TIFFGetField (tiff_handle, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
	TIFFGetField (tiff_handle, TIFFTAG_PHOTOMETRIC, &photometric);
	TIFFGetField (tiff_handle, TIFFTAG_SAMPLEFORMAT, &sampleformat);

    tiff->xmax = xsize;
    tiff->ymax = ysize;

    tiff->data = (float**)calloc(tiff->xmax * tiff-> ymax, 4);

    TIFFGetField(tiff_handle, TIFFTAG_PLANARCONFIG, &planarconfig);
    TIFFReadContigStripData(tiff_handle,4,bitspersample,tiff->data);
    TIFFClose(tiff_handle);
}

void write_tiff(const char *filename, img_t *tiff) {
}