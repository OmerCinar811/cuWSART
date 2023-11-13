#include "wsart_support.h"

int TIFFReadContigStripData(TIFF* tif, int elemsize, int bpp, void* dataptr) {
    unsigned char *buf;
    tsize_t scanlinesize = TIFFScanlineSize(tif);
    long tiffstripsize;
    unsigned char* p;
    unsigned short* ps;
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
			else if ((bpp==16) && (elemsize==6))	/* 48-bit color RGB */
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
    uint16_t bitspersample, samplesperpixel, planarconfig, photometric, sampleformat;
    
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

    tiff->data = (float**)calloc(tiff->xmax * tiff->ymax, 4);


    TIFFGetField(tiff_handle, TIFFTAG_PLANARCONFIG, &planarconfig);
    // printf("\nplanarconfig = %d\n", planarconfig);
    TIFFReadContigStripData(tiff_handle,4,bitspersample,tiff->data);
    TIFFClose(tiff_handle);
}


int write_tiff(const char *fname, img_t *img) {
    
    TIFF *tif;
    char* p;

	tif = TIFFOpen(fname, "w");
	if (!tif) {
        printf("\nCould not open tif during write\n");
        return -1;
    }
	/* Let's start with some general tags */

	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, img->xmax);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, img->ymax);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, 0);

	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, (int)2);
	TIFFSetField(tif, TIFFTAG_XRESOLUTION, 1200);
	TIFFSetField(tif, TIFFTAG_YRESOLUTION, 1200);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,  32);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,    PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,   SAMPLEFORMAT_IEEEFP);
    
    p = (char*)img->data;
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, img->ymax);
	TIFFWriteRawStrip(tif, 0, p, 4*img->xmax*img->ymax);

	TIFFClose(tif);

	return 0;
}

void copyimg(img_t *in, img_t *out) {
    out->xmax = in->xmax;
    out->ymax = in->ymax;
    out->data = (float**)calloc(out->xmax*out->ymax, 4);
    memcpy(out->data, in->data, 4*out->xmax*out->ymax);
}


/* Methods to make sinogram */
// Assume parallel beam with 1 degree ingrements //
// default interporlate is bilinear
// z is 1
// arraywidth is 512 for    detector arr
void make_sinogram(img_t *img){
    float deltad;
    float deltaphi, phi;
    img_t sinogram;
    float s1x,s1y,s2x,s2y;		/* Source and detector midpoints; offset */
    int a,d,t,line,xm,ym,nv,aw2;
    float* p;
    double buf, atten;
    double cphi, sphi;

	/* Get memory for the sinogram image */
	
	deltaphi = 1;					/* Angular increment in degrees */
	nv = 360.0;	/* Number of views */
	xm = img->xmax; ym = img->ymax;		/* we use that all over the place */
	aw2 = ARRWID/2;					/* Midpoint of detector array */
	d = MIN(xm,ym);						/* Diameter of largest circle completely inside box */
	deltad = (float)d/(float)ARRWID;	/* Size of one detector element in image pixels */

	allocate_image (&sinogram, ARRWID, nv, FLOATT);
	p = (float*)sinogram.data;
	
    line=0;
	for (phi=0; phi<360; phi+=deltaphi)	/* One full revolution including redundancy */
	{
		/* We can pre-compute cos phi and sin phi for this angle */

		cphi = COS(phi);
		sphi = SIN(phi);
		/* Now compute one projection (one view) */

		for (a=0; a<ARRWID; a++)		/* Run along the detector array */
		{
			atten=0;						/* Cumulative attenuation along ray */
			for (t=0; t<d; t++)				/* Run along one line */
			{
				/* convert a and t into x,y of an unrotated image coordinate system */

				s1x = t-d/2;				/* horizontal offset from image center */
				s1y = (a-aw2)*deltad;		/* Vertical offset of ray from image center */

				/* Rotate s1x,s1y into the image coordinate system */

				s2x = xm/2 + s1x*sphi - s1y*cphi;
				s2y = ym/2 -s1x*cphi - s1y*sphi;

				/* add image value at that point to total attenuation */

				if ((s2x>=0) && (s2x<xm) && (s2y>=0) && (s2y<ym))
				{
					buf = ireadbuf (img, s2x, s2y);
					atten += buf;
				}
                //printf("phi is %f, a is %d, t is %d, atten is %f\n", phi, a, t, atten);

			}
            //printf("phi is %f, a is %d, t is %d, atten is %f\n", phi, a, t, atten);

			/* Place attenuation in sinogram buffer */
            p[a + (line*sinogram.xmax)] = atten;
		}
		line++;
	}

	handback (img, &sinogram);
	//if (!macrorun) reset_progress ();

} // make_sinogram

void allocate_image(img_t *img, int x, int y, int type) {
    img->xmax = x;
    img->ymax = y;
    printf("got here and x and y is %d %d\n", x, y);
    if(type == 1) { //float
        img->data = (float**)calloc(x*y, 4);    
    } else if(type == 2) { //rgb triplet
        //img->data = (float**)calloc(x*y, 3);
    }
    //img->data = (float**)calloc(x*y, 4);

    if(img->data == NULL) {
        printf("Could not allocate image\n");
        exit(1);
    }
} // allocate_image

double ireadbuf(img_t *img, double x, double y) {
    //double a1,a2,a3,a4,a5,a6;	/* 4 nearest neighbors */
    int x1;
    int y1;

    // x1=floor(x); x2=floor(x+1);	/* Determine the 4 neighbors */
    // y1=floor(y); y2=floor(y+1);
    // a1 = readbuf_flt (img,x1,y1);
    // a2 = readbuf_flt (img,x2,y1);
    // a3 = readbuf_flt (img,x1,y2);
    // a4 = readbuf_flt (img,x2,y2);
    x1 = floor (0.5+x);
    y1 = floor (0.5+y);
    printf("x1 is %d and y1 is %d\n", x1, y1);
    // exit(0);
    return readbuf_flt (img, x1,y1);
    
    // x = x-x1; y = y-y1;	/* Should always be in range 0...1 */
    // a5 = (1-x)*a1 + x*a2;
    // a6 = (1-x)*a3 + x*a4;
    // return (1-y)*a5 + y*a6;
} // ireadbuf

double readbuf_flt (img_t *img, int x, int y) {
    float* p4;

	if (x<0) x=0; else if (x>=img->xmax) x=img->xmax-1;
	if (y<0) y=0; else if (y>=img->ymax) y=img->ymax-1;
 
	p4=(float*)img->data;
    float temp = p4[x+img->xmax*(y+img->ymax)];
    printf("Temp is %f\n", temp);
    return p4[x+img->xmax*(y+img->ymax)];
	
} // readbuf_flt

void handback(img_t *dest, img_t *src) {
	freebuf (dest);
	dest->data = src->data; 
	dest->xmax = src->xmax;
	dest->ymax = src->ymax;
	//make_minmax (dest);
}

void freebuf(img_t *img) {
    if (img->data) free(img->data);
	img->data = NULL;
	img->xmax = img->ymax = 0;
}