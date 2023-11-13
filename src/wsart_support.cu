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
// arraywidth is 512 for detector arr
void make_sinogram(img_t *img){
    float d, deltad;
    float deltaphi, phi;
    img_t sinogram, tracegram;
    float s1x, s1y, s2x, s2y, sdx, sdy;
    int a,cnt,line,column,xm,ym,ix,iy,sc,nv,aw2;
    uint64_t idx;
    float x1,y1,x2,y2,h,dx,dy,x,y;
    float* p;   
    rgbtriplet* cp; 
    rgbtriplet r1,r2;   
    double buf, atten;  
    double *sintbl, *costbl;    

    deltaphi = 1;
    nv = (int) (0.5+360.0/deltaphi);
    xm = img->xmax; ym = img->ymax;
    aw2 = 256; /* 512/2 */

    allocate_image(&sinogram, ARRWID, nv, FLOATT);
    p = (float*)sinogram.data;
    allocate_image(&tracegram, 2*xm, 2*ym, RGBTRPT);
    cp = (rgbtriplet*)tracegram.data;
    d = 0.5*sqrt( SQR(xm) + SQR(ym));

    deltad = 2*d/ARRWID;
    sintbl = (double*)calloc(nv ,sizeof(double));
    costbl = (double*)calloc(nv ,sizeof(double));
    sc = 0;
    for(phi = 0; phi < 360; phi += deltaphi) { // one revolution
        sintbl[sc] = cos (PHI);
		costbl[sc] = -sin (PHI);
		sc++;
    }

    line=0; idx=0; sc=0;
	for (phi=0; phi<360; phi+=deltaphi)	/* One full revolution including redundancy */
	{
/*		s1x = - (d+2)*costbl[sc]; s1y = - (d+2)*sintbl[sc];	   Source midpoint */
		s1x = - d*costbl[sc]; s1y = - d*sintbl[sc];			/* Source midpoint */
		s2x = -s1x; s2y = -s1y;								/* Detector midpoint */
		column=0;

		for (a = 0; a<ARRWID; a++)				/* One detector line */
		{
			sdx = -deltad*((a-aw2)*sintbl[sc]);
			sdy =  deltad*((a-aw2)*costbl[sc]);			/* Detector element offset */
			x2=s2x+sdx; y2=s2y+sdy;					/* Endpoint of beam */
			x1=s1x+sdx; y1=s1y+sdy;				/* Startpoint at source */
			dx = x2-x1; dy= y2-y1;
			h = sqrt (SQR(dx)+SQR(dy));				/* Path length */
			dx /= h; dy /= h;						/* steps for unit step length */
			
			/* Walk along the beam path */
			
			x=x1; y=y1; cnt=0; atten=0.0;
			do
			{
				buf = ireadbuf(img, x+0.5*xm,y+0.5*ym);
				atten += buf;
				
				/* This part generates the X-ray traces in the secondary wnd */
				if ((x+xm>=0) && (x+xm<2*xm) && (y+ym>=0) && (y+ym<2*ym))
				{
					r2 = (rgbtriplet){0,0,0};
					r2.blue = (int) (10.0+240*phi/360);
					r2.red = 255-r2.blue;
					r1=r2; r1.green=255;
					
					ix = (int) (0.5+x); iy = (int) (0.5+y);

					if ((x+0.5*xm>=0) && (x+0.5*xm<xm) && (y+0.5*ym>=0) && (y+0.5*ym<ym))
						cp[ ( (ix+xm) + 2*xm*(iy+ym) )] = r1;
					else
						cp[ ( (ix+xm) + 2*xm*(iy+ym) )] = r2;
				}

				x+=dx; y+=dy; cnt++;
			}
			while ( (SQR(x-x2)+SQR(y-y2)) > 1.0);	/* Proximity of endpoint */
			p[column + (ARRWID)*line] = atten;
			column++;
		}
		line++;
		sc++;
	}


} // make_sinogram

void allocate_image(img_t *img, int x, int y, int type) {
    img->xmax = x;
    img->ymax = y;

    if(type == 1) { //float
        img->data = (float**)calloc(x*y, 4);    
    } else if(type == 2) { //rgb triplet
        img->data = (float**)calloc(x*y, 3);
    }

    if(img->data == NULL) {
        printf("Could not allocate image\n");
        exit(1);
    }
} // allocate_image

double ireadbuf(img_t *img, double x, double y) {
    double a1,a2,a3,a4,a5,a6;	/* 4 nearest neighbors */
    double xval[4];				/* 4 row interpolations in bi-cubic */
    double A[4];				/* Cubic interpolation coefficents */
    int x1,x2,y1,y2,i;

    x1=floor (x); x2=floor (x+1);	/* Determine the 4 neighbors */
    y1=floor (y); y2=floor (y+1);
    a1 = readbuf_flt (img,x1,y1);
    a2 = readbuf_flt (img,x2,y1);
    a3 = readbuf_flt (img,x1,y2);
    a4 = readbuf_flt (img,x2,y2);
    
    x = x-x1; y = y-y1;	/* Should always be in range 0...1 */
    a5 = (1-x)*a1 + x*a2;
    a6 = (1-x)*a3 + x*a4;
    return (1-y)*a5 + y*a6;
} // ireadbuf

double readbuf_flt (img_t *img, int x, int y) {
    float* p4;

	if (x<0) x=0; else if (x>=img->xmax) x=img->xmax-1;
	if (y<0) y=0; else if (y>=img->ymax) y=img->ymax-1;
 
	p4=(float*)img->data;
    return p4[x+img->xmax*(y+img->ymax)];
	
} // readbuf_flt