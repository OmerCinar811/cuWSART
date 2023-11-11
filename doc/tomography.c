
#include <gtk/gtk.h>
#include <gdk/gdk.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <fftw3.h>

#include "includeall.h"



/* Redefine fftw 2's fftw_real for fftw3: */

#ifndef fftw_real

#define fftw_real 	double

#endif

/* In addition, the struct containing re and im no longer exists. Use
   these macros instead: */

#define c_re(C)	((C)[0])
#define c_im(C)	((C)[1])





/**************************************************************************

 #   #     #     #   #   #####  
 ## ##    # #    #  #    #      
 # # #   #   #   # #     #      
 # # #   #   #   ##      ####   
 #   #   #####   # #     #      
 #   #   #   #   #  #    #      
 #   #   #   #   #   #   #####  

 
  ####     *#*    #   #    ####     ####    ####      #     #   #  
 #    #     #     #   #   #    #   #    #   #   #    # #    ## ##  
 #          #     ##  #   #    #   #        #   #   #   #   # # #  
  ####      #     # # #   #    #   #   ##   ####    #   #   # # #  
      #     #     #  ##   #    #   #    #   # #     #####   #   #  
 #    #     #     #   #   #    #   #    #   #  #    #   #   #   #  
  ####     *#*    #   #    ####     ####    #   #   #   #   #   #  


****************************************************************************/

/*	Functions that simulate the forward projection.

	make_sinogram					Create a sinogram through ray tracing
	make_sinogram_simple			Faster version of the ray tracing algorithm, some fixed assumptions
	make_sinogram_radon				Create a sinogram through Fourier slice theorem (defunct?)




/* This little fun function makes a sinogram of the source
   image and returns it in img. We make the following assumptions:
   
   - The image is xm wide and ym high
   - We use MIN(x,y) as the length of the detector array
   - The midpoint of the detector array is at aw2=arraywidth/2
   - The center of the image coordinate system is at xm/2, ym/2
   - Rotation occurs around the angle theta
   
   - Let d be the half-diagonal, d=sqrt(sqr(x)+sqr(y))
   - Let S1 be the midpoint of the source and S2 the midpoint of the detector
   - S2 has the distance d from the center with S2=(x2,y2) and
     S2= (d*cos theta, d*sin theta)
   - S1 is in the opposite direction, thus S1= -S2
   - Each individual sensor element is of thickness delta d = d/256 (512-array)
   - The k-th element is offset from S2 by
     delta x2 = k * delta d * sin theta - delta d / 2
	 delta y2 = k * delta d * cos theta - delta d / 2
	 (k runs from -255 to +256 for a total of 512 elements)
   - In parallel beam geometry, S1 is offset by the same amount
   - In fan beam geometry, S1 is constant
   - Everything outside of the image rectangle has a value of zero
   
*/

/* phi in radians */
#define PHI		(phi*0.017453293)


void make_sinogram (image_type* img, int z, int arraywidth, float anginc, int geometry)
{
float d, deltad;
float deltaphi, phi;
image_type sinogram, tracegram;
float s1x,s1y,s2x,s2y,sdx,sdy;		/* Source and detector midpoints; offset */
int a,cnt,line,column,xm,ym,ix,iy,sc,nv,aw2;
long idx;
float x1,y1,x2,y2,h,dx,dy,x,y;
float* p;
rgbtriplet* cp;
rgbtriplet r1,r2;
double buf, atten;
double *sintbl, *costbl;


	/* Get memory for the sinogram image */
	
	deltaphi = anginc;
	nv = (int) (0.5+360.0/deltaphi);	/* Number of views */
	xm = img->xmax; ym = img->ymax;		/* we use that all over the place */
	aw2 = arraywidth/2;					/* Midpoint of detector array */

	if (allocate_image (&sinogram ,B_FLOAT, arraywidth, nv, 1)==-1)
		fatal_memerr();
	p = sinogram.data;

	if (allocate_image (&tracegram ,B_RGB, 2*xm, 2*ym, 1)==-1)
		fatal_memerr();
	cp = tracegram.data;


	/* Choose one of the following two. The first one takes the maximum
	   diagonal dimensions of the image to determine the reconstruction
	   circle radius with the following properties:
	   	- includes _all_ of the image
		- shrinks image size by 71% after reconstruction
	   The second one takes the smaller box dimension as reconstruction
	   circle diameter with the following properties:
	   	- does not consider parts of the image outside the circle (i.e. corners)
		- image size is preserved after reconstruction
	*/

	d = 0.5*sqrt( SQR(xm) + SQR(ym));	/* half the diagonal width */

	/* or */
	
	d = 0.5*MIN(xm,ym);					/* largest circle completely inside box */

	deltad = 2*d/arraywidth;	/* Size of one detector element */

	/* For faster computation, we prepare sin/cos tables in advance.
		In the inner loop, phi is i*deltaphi, so each table element
		at [i] contains the value sin(phi*deltaphi)
		Note - for compatibility with ART projection, we rotate the
		image by 90 degrees ccw - that's why sintbl and costbl are 
		not exactly what their name implies... */
	
	sintbl = calloc (nv, sizeof(double));
	costbl = calloc (nv, sizeof(double));
	sc=0;
	for (phi=0; phi<360; phi+=deltaphi)	/* One full revolution including redundancy */
	{
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

		for (a = 0; a<arraywidth; a++)				/* One detector line */
		{
			sdx = -deltad*((a-aw2)*sintbl[sc]);
			sdy =  deltad*((a-aw2)*costbl[sc]);			/* Detector element offset */
			x2=s2x+sdx; y2=s2y+sdy;					/* Endpoint of beam */
			if (geometry==0)						/* Parallel-beam */
			{
				x1=s1x+sdx; y1=s1y+sdy;				/* Startpoint at source */
			}
			else									/* fan-beam */
			{
				x1=s1x; y1=s1y;						/* Startpoint at source */
			}
			dx = x2-x1; dy= y2-y1;
			h = sqrt (SQR(dx)+SQR(dy));				/* Path length */
			dx /= h; dy /= h;						/* steps for unit step length */
			
			/* Walk along the beam path */
			
			x=x1; y=y1; cnt=0; atten=0.0;
			do
			{
				buf = ireadbuf (img, x+0.5*xm,y+0.5*ym,z,intp_pref);
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
			p[column + (arraywidth)*line] = atten;
			column++;
		}
		line++;
		sc++;
	}

	handback (img, &sinogram);

	make_minmax (&tracegram);
	secondary_display (&tracegram,0,"Ray traces");

	freebuf (&tracegram);
	free (sintbl); free (costbl);

}


/* A variation of the above. By not considering geometric expansion, we can simplify
	the sinogram computation dramatically. Consider a projection along the x axis.
	In that case, the detector position (vertical) and the line trace
	(horizontal) are Cartesian coordinates in the non-rotated image. We
	simply rotate these coordinates to get the image coordinates of the 
	rotated image (with an inverse rotation matrix to reflect the rotation of
	the projection wrt image) */


void make_sinogram_simple (image_type* img, int z, int arraywidth, float anginc, int fanbeam_dummy)
{
float deltad;
float deltaphi, phi;
image_type sinogram;
float s1x,s1y,s2x,s2y;		/* Source and detector midpoints; offset */
int a,d,t,line,xm,ym,sc,nv,aw2;
float* p;
double buf, atten;
double cphi, sphi;
long idx;


	/* Get memory for the sinogram image */
	
	deltaphi = anginc;					/* Angular increment in degrees */
	nv = (int) (0.5+360.0/deltaphi);	/* Number of views */
	xm = img->xmax; ym = img->ymax;		/* we use that all over the place */
	aw2 = arraywidth/2;					/* Midpoint of detector array */
	d = MIN(xm,ym);						/* Diameter of largest circle completely inside box */
	deltad = (float)d/(float)arraywidth;	/* Size of one detector element in image pixels */

	if (allocate_image (&sinogram ,B_FLOAT, arraywidth, nv, 1)==-1)
		fatal_memerr();
	p = sinogram.data;

	line=0; idx=0; sc=0;
	for (phi=0; phi<360; phi+=deltaphi)	/* One full revolution including redundancy */
	{
		/* We can pre-compute cos phi and sin phi for this angle */

		cphi = COS(phi);
		sphi = SIN(phi);

		dp (1, "Make sinogram: Processing projection at angle %5.1f\n",phi);
		if (!macrorun) indicate_progress (phi/360.0);

		/* Now compute one projection (one view) */

		for (a=0; a<arraywidth; a++)		/* Run along the detector array */
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
					buf = ireadbuf (img, s2x, s2y, z, intp_pref);
					atten += buf;
				}

			}

			/* Place attenuation in sinogram buffer */

			p[a + line*sinogram.xmax] = atten;
		}
		line++;
	}

	handback (img, &sinogram);
	if (!macrorun) reset_progress ();

}





/* Alternatively, we may use the Fourier slice theorem for sinogram generation.

	- Transform the input image
	- Grab lines from the transform subtending the detector angle theta
	- Arrange them in a sinogram-like array
	- Inverse-transform that array line-by-line
	
*/


void make_sinogram_radon (image_type* img, int z, int arraywidth, float anginc)
{
image_type sino;		/* the final sinogram */
image_type dbgimg;
int x,y, xm, ym, nv, deltaphi,i,j, view;
int xh, yh, exts;
float *p; float *dbgp;
float deltax, deltay, ry, rx, theta;
double tempr, tempi;
fftw_complex *xform;		/* Primary 2D FFT of original image */
fftw_complex *lp;			/* One detector line */
fftw_complex *lpfd;			/* One detector line, freq domain */
fftw_complex *in;			/* 2D cplx-valued input to the DFT */
fftw_plan plan1d, plan2d;


	/* If the number of detectors is odd, make it even */
	
	if (arraywidth & 1) arraywidth+=1;

	deltaphi = anginc;
	nv = (int) (0.5+360.0/deltaphi);	/* Number of views */
	xm = img->xmax; ym = img->ymax;	/* we use that all over the place */
	xh=xm/2; yh=ym/2;
	exts = MIN(xm,ym);				/* diameter of the reconstruction circle */

	if (allocate_image (&sino ,B_FLOAT, arraywidth, nv, 2)==-1)
		fatal_memerr();
	p = sino.data;

	if (allocate_image (&dbgimg ,B_FLOAT, xm,ym, 2)==-1)
		fatal_memerr();
	dbgp = dbgimg.data;

	in = (fftw_complex*) fftw_malloc (2*xm*ym*sizeof (fftw_complex));
	if (!in) fatal_memerr();

	xform = (fftw_complex*) fftw_malloc (2*xm*ym*sizeof (fftw_complex));
	if (!xform) fatal_memerr();

	lp = (fftw_complex*) fftw_malloc ((arraywidth+1)*sizeof (fftw_complex));
	if (!lp) fatal_memerr();

	lpfd = (fftw_complex*) fftw_malloc ((arraywidth+1)*sizeof (fftw_complex));
	if (!lpfd) fatal_memerr();


	plan1d = fftw_plan_dft_1d (arraywidth, lpfd, lp, FFTW_BACKWARD, FFTW_ESTIMATE);
	plan2d = fftw_plan_dft_2d (xm,ym, in, xform, FFTW_FORWARD, FFTW_ESTIMATE);


	/* Let's 2D-transform the original image */

	dp (2, "Sinogram by Radon: Forward FFT...\n");	

	for (y=0; y<xm; y++)
		for (x=0; x<xm; x++)
		{
			c_re(in[x*ym+y]) = readbuf_flt (img,x,y,z);
			c_im(in[x*ym+y]) = 0;
		}

	for (y=0; y<xm; y++)
		for (x=0; x<xm; x++)
		{
			c_re(xform[(y + ym*x)]) = 0.0;	/* Real component */
			c_im(xform[(y + ym*x)]) = 0.0;	/* Imaginary component */
		}
			
	
	fftw_execute (plan2d);
	dp (2, "...done.\n");
	
	/* Swap quadrants. We can reuse the "in" buffer.
		Visualize the 2D transform. Note that after dft_r2c, only the 
	   	y columns from 0 to yh+1 are filled. Note also that the output
	   	is scaled, so we need to unscale (we do this after the forward FT here) */
	
	for (y=0; y</*(yh+1)*/ym; y++)
		for (x=0; x<xm; x++)
		{

			if (x<xh) i=x+xh; else i=x-xh;
			if (y<yh) j=y+yh; else j=y-yh;
			tempr = c_re(xform[y+ym*x])/(xm*ym);	/* Scaled real component */
			tempi = c_im(xform[y+ym*x])/(xm*ym);	/* Scaled imaginary component */

			dbgp[i+xm*j] = tempr;
			dbgp[i+xm*j+xm*ym] = tempi;
			c_re(in[j+ym*i]) = tempr;			/* Swapped; DC in the center */
			c_im(in[j+ym*i]) = tempi;
		}
	make_minmax (&dbgimg);
	


	/* View by view, grab one line from the 2D transform,
	   inverse transform, and enter into sino */

	dp (2, "Sinogram by Radon: Inverse FFT...\n");	

	view=0;
	for (theta=0; theta<180; theta+=anginc)
	{
		deltax = exts*cos(theta*PI/180)/arraywidth;
		deltay = exts*sin(theta*PI/180)/arraywidth;
		rx=0; ry=0;

		/* Positive frequencies starting with DC at lp[0] */

		for (i=0; i<arraywidth/2; i++)
		{
			x = (int) (0.4999+xh+rx); y = (int)(0.4999+yh+ry);
			if (i>0)
			{
				c_re(lpfd[arraywidth-i]) = c_re(in[x*ym+y]);
				c_im(lpfd[arraywidth-i]) = c_im(in[x*ym+y]);
			}
			c_re(lpfd[i]) = c_re(in[x*ym+y]);
			c_im(lpfd[i]) = c_im(in[x*ym+y]);

			x = (int) (0.4999+xh-rx); y = (int)(0.4999+yh-ry);
			c_re(lpfd[arraywidth-i-1]) = c_re(in[x*ym+y]);
			c_im(lpfd[arraywidth-i-1]) = c_im(in[x*ym+y]);

			rx += deltax;
			ry += deltay;
		}

		/* inverse-transform lp */
		
		fftw_execute (plan1d);
		
		/* Enter the result into the sino image */
		
		for (i=0; i<arraywidth; i++)
		{
			p[i + view*arraywidth] = c_re(lpfd[i]);
			p[i + view*arraywidth + arraywidth*nv] = c_im(lpfd[i]);

		}
		view++;
	
	}

	make_minmax (&dbgimg);
	vsrtwrite (&dbgimg, "radonfft.vflt",0,1);
	secondary_display (&dbgimg, 0, "2D FFT");



	dp (2, "...done.\n");


	/* Finish off: Free memory, hand back imge */

	fftw_free (in);
	fftw_free (xform);
	fftw_free (lp);
	fftw_free (lpfd);
	fftw_destroy_plan (plan2d);
	fftw_destroy_plan (plan1d);

	make_minmax (&sino);

	freebuf (&dbgimg);

	handback (img, &sino);

}


/* Here, we simulate a PET sinogram by Monte Carlo simulation. We assume two "hot spots"
	at (303,165) and at (237,338) in a 512x512 image. Their emission is approximately 
	Gaussian, and we choose for each emission event a random angle for the line of
	coincidence. The distance from the origin is given. Each event gets recorded
	as (distance,angle). We have approximately homogeneous background radiation */

/* Presently, img is not used. We can later refine this by using img as absorbance map */


void make_pet_sinogram (image_type* img, int z, int arraywidth, float anginc)
{
float d, prob, rho;
float deltaphi, phi;
image_type sinogram;
float s1x,s1y,s2x,s2y,s1sd,s2sd;		/* Source and detector midpoints; stdev */
int cnt,xm,ym,sc,nv,aw2,det;
float x,y;
float* p;


	/* Get memory for the sinogram image */
	
	deltaphi = anginc;
	nv = (int) (0.5+180.0/deltaphi);	/* Number of views */
	xm = arraywidth; ym = arraywidth;		/* we use that all over the place */
	aw2 = arraywidth/2;					/* Midpoint of detector array */

	if (allocate_image (&sinogram ,B_FLOAT, arraywidth, nv, 1)==-1)
		fatal_memerr();
	p = sinogram.data;

	d = aw2;					/* largest circle completely inside box */

	/* Define the hot-spot centers scaled to image size */

	s1x = (303.0-256.0)*xm/512.0;			/* x distance from origin */
	s1y = (165.0-256.0)*ym/512.0;			/* y distance from origin */
	s2x = (237.0-256.0)*xm/512.0;			/* x distance from origin */
	s2y = (338.0-256.0)*ym/512.0;			/* y distance from origin */
	s1sd = 5.0*d/256.0;					/* standard deviation wrt 512x512 image */
	s2sd = 10.0*d/256.0;

	cnt=0;
	while (cnt < 50000)				/* Number of events total */
	{
		/* First, we decide whether we have background or hot spot */

		prob = nr_rand2();			/* Uniform 0...1 distribution */

		if (prob < 0.01)				/* Choose background event */
		{
			phi = nr_rand2 ()*180.0;		/* Random angle between 0 and 180 */
			x = gauss_rnd ()*2*d;			/* Random distance from origin, kinda Gaussian */
			y = gauss_rnd ()*2*d;
		}
		else if (prob < 0.3)		/* Weaker upper hot spot */
		{
			phi = nr_rand2 ()*180.0;		/* Random angle between 0 and 180 */
			x = s1x + gauss_rnd()*s1sd;
			y = s1y + gauss_rnd()*s1sd;
		}
		else						/* Stronger lower hot spot */
		{
			phi = nr_rand2 ()*180.0;		/* Random angle between 0 and 180 */
			x = s2x + gauss_rnd()*s2sd*0.8;
			y = s2y + gauss_rnd()*s2sd*1.5;	/* This spot is elliptical */
		}

		/* Here, we have originating point and angle of LOC. Compute rho,
			then one-up the sinogram. */

		rho = x*COS(phi) + y*SIN(phi);
		det = (int)(0.4999 + rho + aw2);
		sc = (int) (0.4999 + phi/deltaphi);
		if ((sc>=0) && (sc < nv) && (det >=0) && (det < arraywidth))
		{
			p[det + sc*arraywidth] += 1.0;
			cnt++;
		}
	}


	handback (img, &sinogram);


}




/************************************************************

	   #      ######    #######
	  # #     #     #      #
	 #   #    #     #      #
	#     #   ######       #
	#######   #   #        #
	#     #   #    #       #
	#     #   #     #      #

************************************************************/

/* Arithmetic reconstruction technique. Projection & reconstruction. */

/*
	OVERVIEW OF ART FUNCTIONS

	find_prime			Find the nearest smaller prime number (needed for projection skew)
	adapt_underrelaxation 	Adaptive underrelaxation (update of underrelaxation coefficient)
	in_recon_circle		returns whether a point (x,y) is inside the reconstruction circle
	hamming				returns the value of a Hamming window at point idx in a total of N points
	Wij_old				Calculate Wij through line intersections, explicit, inefficient, old algorithm
	Wij1				Calculate Wij through line intersections
	precompute_art_area	Helper function for Wij2
	Wij2				Calculate Wij through finite beam/circle overlap area
	Wij_sched			Selector function for Wij1 or Wij2
	art_tricks			Adaptive smoothing trick and constraint trick
	refraction_trick	Interpolate over photon-starved regions of a projection
	art_project			Forward projection: Compute sinogram from image using the Wij
	art_recon			Core of the SART and ART algorithms
	wijtest				A test function fot eh Wij computation. No longer sure how it works...
*/




/* Compute the Wij weight parameters.
	x and y are the coordinates of pixel i
	j and theta are the projection coordinates
	N is the number of detectors AND the image size (NxN)

	We have two versions of Wij. First, we assume a line beam, and we check the
	intersection of each (square) pixel with the line. The length of the line segment
	intersecting the square is our Wij for pixel i and projection j.

	The second version is (hopefully) more accurate as it uses the beam
	overlap with a circular approximation of the pixel. In fact, with a
	suitable pre-computed table, we can even consider Gaussian beam profiles.

	One way for speeding things up is to use implicit multiplications
	for determinant computation rather than the function call. We could
	make use of the constant zeroes to eliminate some multiplications.
*/

/* Choose float or double for ART. NOTE: Float is not sufficiently precise
	to calculate the Wij. With float, we frequently get single or triple
	intersection points. Double avoids this. */

#define AFLOAT	double


int art_use_smtrick, art_use_ctrick, art_use_hamming, art_use_refcorr;
int art_maxiter; float art_epsmax;

AFLOAT *art_sintbl, *art_costbl;			/* Pre-computed sin(theta) and cos(theta) */
AFLOAT *art_overlap;						/* Precomputed circle/beam overlap */



/* Find the nearest prime number that is <= n  */

int find_prime (int n)
{
int primes[30] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
					71, 73, 79, 83, 89, 97, 101, 103, 107, 109};
int i;

	i=29;
	while ((i>0) && (n<primes[i])) i--;
	return MAX (1, primes[i]);
}


/* Adaptive underrelaxation, adjusted after every iteration.
	Underrelaxation coefficient is increased at least 0.01 every time.
	The increase becomes stronger as the normalized error decreases.
	Underrelaxation is clamped to 1.00 (no underrelaxation) */

AFLOAT	adapt_underrelaxation (AFLOAT lambda, AFLOAT epsilon)
{
	if (lambda < 1.0)
	{
		if (epsilon < 0.5) lambda=1.0;
		else if (epsilon < 2.0) lambda += 0.1;
		else if (epsilon < 5.0) lambda += 0.05;
		else lambda += 0.01;
	}
	lambda = MIN(1.0, lambda);
	return lambda;
}



/* Test whether a coordinate x,y is within the reconstruction circle
	of diameter d */


unsigned char in_recon_circle (int x, int y, int d)
{
	d = (d*d) >> 2;		/* now radius square */
	if ((x*x+y*y) > d) return 0;
	return 1;
}

/* Hamming or Raised Cosine window */

AFLOAT hamming (float idx, float N)
{
	return 0.54 + 0.46*cos(TWOPI*idx/N);
}


/* This first function is the explicit one that uses determinants
	to compute the line intersection with the pixel square.
	-- use for testing only. */

MFLOAT Wij_old (int x, int y, int j, float theta, int N)
{
MFLOAT s,deta;
MFLOAT A[4], B[4], C[4];	/* Matrices: A denominator; B and C enumerators for intersection pt */
MFLOAT intx[4], inty[4];	/* Up to (theoretically) four intersection points */
int intpts;					/* no. of intersection pts - should always be either 0 or 2 */

	/* Step 1. The detector distance from the detector array center
		(identical to the ray distance from the origin).
		s and theta define Hesse's normal form of the ray */
	
	if (j >= N/2) 
		s = - (0.5+j-N/2);
	else
		s = 0.5+N/2-1-j;

	/* Step 2: Check four possible intersection points with the
		walls of the pixel.
		Step 2a: Left and right delimiters */

	intpts=0;

	A[0]=1.0; A[1]=0.0;
	A[2]=cos(theta); A[3] = sin(theta);
	deta = det2x2(A);
	if (deta!=0.0)			/* If this det is zero, the rays are parallel */
	{
		/* Compute the intersection points */
		
		B[0]=0; B[1]=N/2-x;
		B[2]=sin(theta); B[3] = -s;
		C[0]=N/2-x; C[1]=1;
		C[2]=-s; C[3] = cos(theta);
		intx[intpts] = det2x2(B)/deta;
		inty[intpts] = det2x2(C)/deta;
		
		/* Are the intersection points inside the pixel boundary? */
	
		if ((inty[intpts] >= (y-N/2)) && (inty[intpts]<=(1+y-N/2)))
			intpts++;

		B[1]=N/2-x-1;
		C[0]=N/2-x-1;
		intx[intpts] = det2x2(B)/deta;
		inty[intpts] = det2x2(C)/deta;
		
		/* Are the intersection points inside the pixel boundary? */
	
		if ((inty[intpts] >= (y-N/2)) && (inty[intpts]<=(1+y-N/2)))
			intpts++;

	}

	/* Step 2b: Top and bottom delimiters */

	A[0]=0.0; A[1]=1.0;
	A[2]=cos(theta); A[3] = sin(theta);
	deta = det2x2(A);
	if (deta!=0.0)			/* If this det is zero, the rays are parallel */
	{
		/* Compute the intersection points */
		
		B[0]=1; B[1]=N/2-y;
		B[2]=sin(theta); B[3] = -s;
		C[0]=N/2-y; C[1]=0;
		C[2]=-s; C[3] = cos(theta);
		intx[intpts] = det2x2(B)/deta;
		inty[intpts] = det2x2(C)/deta;
		
		/* Are the intersection points inside the pixel boundary? */
	
		if ((intx[intpts] >= (x-N/2)) && (intx[intpts]<=(1+x-N/2)))
			intpts++;

		B[1]=N/2-y-1;
		C[0]=N/2-y-1;
		intx[intpts] = det2x2(B)/deta;
		inty[intpts] = det2x2(C)/deta;
		
		/* Are the intersection points inside the pixel boundary? */
	
		if ((intx[intpts] >= (x-N/2)) && (intx[intpts]<=(1+x-N/2)))
			intpts++;
	}

	if (intpts==0) return 0.0;		/* No intersection -> wij=0 */		

	if ((intpts==1) || (intpts==3))
		dp (1, "Error: %d intersection points at (%d,%d), j=%d, theta=%f\n",intpts,x,y,j,theta);
	
	return  sqrt(ARITHM (intx[1]-intx[0], inty[1]-inty[0]));

}


/* This newer algorithm makes use of sin/cos tables and is generally faster */



AFLOAT Wij1 (int x, int y, int j, int thetaidx, int N)
{
AFLOAT intx[4], inty[4];	/* Up to (theoretically) four intersection points */
int intpts;					/* no. of intersection pts - should always be either 0 or 2 */
AFLOAT s,d1,d2,d3;
int xi,yi;					/* shifted coord origin */

	/* Step 1. The detector distance from the detector array center
		(identical to the ray distance from the origin).
		s and theta define Hesse's normal form of the ray */
	
	if (j >= N/2) 
		s = - (0.5+j-N/2);
	else
		s = 0.5+N/2-1-j;
	
	xi = x-N/2;
	yi = y-N/2;

	/* Step 2: Check four possible intersection points with the
		walls of the pixel.
		Step 2a: Left and right delimiters */

	intpts=0;
	if (!F_ISZERO(art_sintbl[thetaidx])) 			/* If this is zero, the rays are parallel */
	{
		/* Compute the intersection points and check if inside pixel boundary */
		
		intx[intpts] = xi;
		inty[intpts] = (s-xi*art_costbl[thetaidx])/art_sintbl[thetaidx];
		if ((inty[intpts] >= yi) && (inty[intpts]<=(1+yi)))
			intpts++;

		intx[intpts] = xi+1;
		inty[intpts] = (s-(xi+1)*art_costbl[thetaidx])/art_sintbl[thetaidx];
		if ((inty[intpts] >= yi) && (inty[intpts]<=(1+yi)))
			intpts++;

	}

	/* If we (coincidentially) found two intersection pts, we can leave 
		and skip the next section for higher efficiency */
	
	if (intpts==2)
		return  sqrt(ARITHM (intx[1]-intx[0], inty[1]-inty[0]));

	/* Step 2b: Top and bottom delimiters */

	if (!F_ISZERO(art_costbl[thetaidx])) 			/* If this is zero, the rays are parallel */
	{
		/* Compute the intersection points and check if inside pixel boundary */
		
		intx[intpts] = (s-yi*art_sintbl[thetaidx])/art_costbl[thetaidx];
		inty[intpts] = yi;
		if ((intx[intpts] >= xi) && (intx[intpts]<=(1+xi)))
			intpts++;
		intx[intpts] = (s-(yi+1)*art_sintbl[thetaidx])/art_costbl[thetaidx];
		inty[intpts] = yi+1;
		if ((intx[intpts] >= xi) && (intx[intpts]<=(1+xi)))
			intpts++;
	}

	if (intpts==0) return 0.0;		/* No intersection -> wij=0 */		

	if ((intpts==1) || (intpts>=3))
		dp (2, "Error: %d intersection points at (%d,%d), j=%d, theta-idx=%d\n",intpts,x,y,j,thetaidx);
	
	if (intpts==1) return 0.0;		/* One intersection - likely it's really zero */

	d1 = ARITHM (intx[1]-intx[0], inty[1]-inty[0]);

	if (intpts==2)
		return  sqrt(d1);
		
	/* This leaves us with the most complicated case. I think we most likely
		hit a corner, so we need to find the maximum of the three distances */

	d2 = ARITHM (intx[2]-intx[0], inty[2]-inty[0]);
	d3 = ARITHM (intx[2]-intx[1], inty[2]-inty[1]);
	
	if ((d1>d2) && (d1>=d3))
		return (sqrt(d1));
	if ((d2>d3) && (d2>=d1))
		return (sqrt(d2));
	return (sqrt(d3));

}


/**************************

	Circular approximation Wij

	Calculate the weight coefficient Wij for a circular approximation of the pixel square.
	The circle has a radius of r, which should usually be 0.5 as we assume unit pixel size.
	The pixel center (and thus the center of the circle) is at (s,t) where s is parallel to
	the beam and t is perpendicular to it. Clearly, s is of no relevance for the Wij.

	The beam has a centerline distance t=rho from the origin (and center of rotation).
	The beam's support (=thickness) is 2u (u to each side of rho). Normally, we use u=0.5
	as the beam thickness determines our pixel size, but this may vary for pseudoresolution.
	The beam may have a Gaussian profile.

	We can test each pixel P=(x,y). We get the coordinates of P=(s,t) in the rotated
	coordinate system by simple coordinate transformation:

		s = x cos(theta)  + y sin(theta)
		t = -x sin(theta) + y cos(theta)

	Note that we only need to compute the t coordinate. From this, the distance of the beam
	center rho to the circle center t is d=|rho-t|; the magnitude is taken due to symmetry.
	The overlap is zeta = u+r-d, and the circle circumference begins at a distance from the
	beam centerline of h=d-r. Note that if h<0, we cross the beam centerline inside the circle.

******************/


/* Precompute the overlap area. 

	1) For now, we assume a uniform beam profile (otherwise, we'd need a 2D array)
	2) We assume a normalized radius (0 < zeta < 1 for a half circle)
		discretized in 10 steps
	3) This step must be done at the same time we compute the sin/cos tables

	u is the beam half-width, usually 0.5; r the radius (usually 0.5) of the circle.
	Both are usually 0.5. The resulting array, art_overlap[], contains the overlap
	area as a function of zeta, that is, the overlap distance (circle rim to beam edge)
*/


int precompute_art_area (float u, float r)
{
int i;
int N=10;			/* Discrete support points */
AFLOAT x,y;			/* Coordinate points along the circumference */
AFLOAT x1,y1;


	art_overlap = calloc (1+2*N, sizeof (AFLOAT));

	art_overlap[0]=0;				/* Just touching the rim of the circle */
	y=(N-1)*r/N; x=sqrt(SQR(r)-SQR(y));		/* First point on circumference */

	/* The first segment is approximated by a triangle, and we have to re-use (x,y) in the next iteration */

	art_overlap[2*N] = 
		art_overlap[1] = x*(r-y);				/* Triangle each side = a square */

	/* The next increments are trapezoidal -- note that we only store the increments
		at this time. We'll integrate further down. */

	for (i=2; i<=N; i++)
	{
		x1=x; y1=y;							/* Save previous point on circumference */
		y -= r/N;							/* Slide down 1/10 of radius */
		x = sqrt(SQR(r)-SQR(y));			/* Corresponding x coord of circumference */
		art_overlap[2*N+1-i] = 
			art_overlap[i] = (x+x1)*r/N;	/* Trapezoid area */
	}

	/* Last, add up the incremental areas to create the overlap area.
		Note the special cases (for higher accuracy) of N (half circle) and 2N */

	for (i=2; i<2*N; i++)
	{
		if (i==N)
			art_overlap[N] = PI*r*r*0.5;
		else
			art_overlap[i] += art_overlap[i-1];
	}

	art_overlap[2*N] = PI*r*r;

	return 2*N;				/* Return the number of data points generated */
}



/* Circular overlap Wij */

AFLOAT Wij2 (int x, int y, int j, int thetaidx, int N)
{
AFLOAT t,d,u,r,val;
AFLOAT rho, zeta;
int xi,yi,i;					/* shifted coord origin */

	/* Step 1. The detector distance from the detector array center
		(identical to the ray distance from the origin).
		rho and theta define Hesse's normal form of the ray.
		Also, compute the center of the reconstruction area (xi,yi) */
	
	if (j >= N/2) 
		rho = - (0.5+j-N/2);
	else
		rho = 0.5+N/2-1-j;
	
	xi = x-N/2;
	yi = y-N/2;
	u = 0.5;				/* Beam half-width, hard-coded for now */
	r = 0.5;				/* Pixel radius, hard-coded for now */

	/* Step 2: Subject the pixel in question to coordinate transform into the (s,t)
		coordinate system. We are only interested in t (perpendicular to beam)
		because it determines the center-center distance of pixel circle and beam. */

	t = yi*art_costbl[thetaidx] - xi*art_sintbl[thetaidx];
	d = fabs(rho-t);										/* Center-center distance */
	zeta = u+r-d;											/* Overlap */
	if ((zeta<0) || (zeta >=2*r))
		return 0;											/* Distance too big, no overlap at all */

	/* Step 3: There is some overlap. Interpolate it from the art_overlap table */

	val = 10.0*zeta/r;
	i = (int)floor(val);
	t = val-i;					/* The fractional value; re-use t, it is no longer needed */
	val = t*art_overlap[i+1] + (1.0-t)*art_overlap[i];
	return val;
}



/* This function is set to schedule either the line-intersection Wij1 or the circle-area Wij2 */

AFLOAT Wij_sched (int x, int y, int j, int thetaidx, int N)
{
	return Wij2 (x,y, j, thetaidx, N);
}



/* *********************************************************

	Tricks are correction steps applied between iterations.
	(See Gabor Herman, "Image Reconstruction from Projections",
	Chapter 11, pg. 193ff, New York: Academic Press, 1980)

	This implementation is a selective smoothing filter (trick==1)
	and a constraint filter (all data must be >=0) for trick==2.
	The image img passed to this function is an iteration of the
	reconstruction space
*/


void art_tricks (image_type *img, char trick)
{
float *p, *q;
int x,y,idx,xm,cnt;
float w,s,t;


	g_assert (img->imgtype==B_FLOAT);
	p = img->data;
	xm = img->xmax;

	if (trick==1)
	{
		/* Create a workspace for the smoothed data */
		q = (float*)calloc (img->xmax*img->ymax, sizeof (float));
 		if (!q) fatal_memerr();

		t = 0.01;						/* Smoothing threshold */
		cnt=0;							/* TODO: Threshold should be adaptive! */

		/* The following is a conditional (adaptive) convolution with the kernel

			|  1  4  1  |
			|  4  9  4  |
			|  1  4  1  |

		   where each neighbor gets included in the convolution *only* if it
		   deviates less than 0.01 from the center pixel. This means that we
		   apply smoothing only to already-smooth regions and this conditional
		   smoothing filter is edge-preserving.
		*/

		for (y=1; y<img->ymax-1; y++)
			for (x=1; x<img->xmax-1; x++)
			{
				idx = x + img->xmax * y;

				w=9.0; s=9.0*p[idx];		/* central pixel */

				/* direct neighbors */

				if ( fabs(p[idx-1]-p[idx]) < t) { w+=4.0; s+=4.0*p[idx-1]; }
				if ( fabs(p[idx+1]-p[idx]) < t) { w+=4.0; s+=4.0*p[idx+1]; }
				if ( fabs(p[idx+xm]-p[idx]) < t) { w+=4.0; s+=4.0*p[idx+xm]; }
				if ( fabs(p[idx-xm]-p[idx]) < t) { w+=4.0; s+=4.0*p[idx-xm]; }

				/* Diagonal neighbors */

				if ( fabs(p[idx-xm-1]-p[idx]) < t) { w+=1.0; s+=p[idx-xm-1]; }
				if ( fabs(p[idx-xm+1]-p[idx]) < t) { w+=1.0; s+=p[idx-xm+1]; }
				if ( fabs(p[idx+xm-1]-p[idx]) < t) { w+=1.0; s+=p[idx+xm-1]; }
				if ( fabs(p[idx+xm+1]-p[idx]) < t) { w+=1.0; s+=p[idx+xm+1]; }

				q[idx] = s/w; if (w>9.0) cnt++;
			}

		for (y=1; y<img->ymax-1; y++)
			for (x=1; x<img->xmax-1; x++)
			{
				idx = x + img->xmax * y;
				p[idx]=q[idx];				/* Copy workspace back into image */
			}
		dp (2, "Smoothing trick: %d of %d pixels corrected\n",cnt,xm*img->ymax);

	}
	else if (trick==2)
	{
		cnt = 0;
		for (idx=0; idx<img->xmax*img->ymax; idx++)
		{
			if (p[idx] < 0.0) { p[idx]=0.0; cnt++; }
		}
		dp (2, "Constraint trick: %d of %d pixels corrected\n",cnt,xm*img->ymax);
	}
}


/* Refraction trick: Interpolate over refraction shadow zones in the sinogram.
	Note that trick==0 does nothing (omission of shadow zones from sinogram) 
	Question -- the algorithm here assumes that there are no negative values
	in the sinogram. Is this acceptable?

	The refraction trick has been published in
	Haidekker MA. Optical Transillumination Tomography with tolerance against 
	refraction mismatch. Computer Programs and Methods in Biomedicine 2005; 80: 225-235. 
*/


void refraction_trick (image_type *img, int trick)
{
float *p, *q;
int theta,t1,t2,t,idx,x,y,cnt;
float buf,buf1,r;


	g_assert (img->imgtype==B_FLOAT);
	p = img->data;

	q = (float*)calloc (img->xmax*img->ymax, sizeof (float));
	if (!q) fatal_memerr();

	if (trick==1)				/* Interpolation trick */
	{
		cnt=0;
		for (theta=0; theta < img->ymax; theta++)	/* Run over all projections */
		{
			t1=0; idx = theta*img->xmax;
			while (t1 < img->xmax-1)
			{
				t2=t1+1;
				buf = p[t1 + idx];
				if (buf >= -1e-8)
				{
					q[t1 + idx] = buf;
					t1++;
				}
				else								/* shadow zone */
				{
					buf = p[t1-1 + idx];
					while ((t2 < img->xmax-1) && (p[t2 + idx] <= -1e-8)) t2++;
					buf1 = p[t2 + idx];
					for (t=t1; t<t2; t++)
					{
						r = (float)(t-t1)/(float)(t2-t1); cnt++;
						q[t + idx] = r*buf1 + (1.0-r)*buf;
					}
					t1=t2;
				}
			}
		}
		for (y=0; y<img->ymax; y++)
			for (x=1; x<img->xmax-1; x++)
			{
				idx = x + img->xmax * y;
				p[idx]=q[idx];
			}
		dp (2, "Refraction trick: %d pixels interpolated.\n",cnt);
	}
	else if (trick==2)				/* Maximum variance trick */
	{
		for (theta=0; theta < img->ymax; theta++)	/* Run over all projections */
		{
			idx = theta*img->xmax;
			for (t=0; t<img->xmax; t++)
			{
				if (p[t1 + idx] < -1e-8)			/* For each shadow zone pixel... */
				{
					r = fabs(img->xmax/2 - t);		/* and....? Now what? */
				}
			}
		}

	}

	free (q);

}





/* Compute projections (sinogram ) using arithmetic projection.
	We could probably increase efficiency incredibly if we swap
	the loops. If x and y are the outer loops, we could cache
	the Wij for one sweep over j,theta */

void art_project (image_type* img, int z, int arraywidth, float anginc)
{
image_type sinogram;
int j,x,y,line,cnt,sc,halfsino, myjob;
AFLOAT sum,nv,buf,phi,deltaphi;
float* p;

#define MIRRORME


	/* Get memory for the sinogram image */
	
	deltaphi = anginc;
	arraywidth &= 0xfffe;				/* Force number of detectors to be even */
	nv = (int) (0.5+360.0/deltaphi);	/* Number of views */
	cnt=0;								/* Silence a compiler warning */

	if (allocate_image (&sinogram ,B_FLOAT, arraywidth, nv, 1)==-1)
		fatal_memerr();
	p = sinogram.data;

	sc=0; halfsino=0; myjob=0;
	if (art_sintbl==NULL)
	{
		art_sintbl = calloc (nv, sizeof(AFLOAT));
		art_costbl = calloc (nv, sizeof(AFLOAT));
		myjob=1;
		for (phi=0; phi<360; phi+=deltaphi)	/* One full revolution including redundancy */
		{
			art_sintbl[sc] = sin (PHI);
			art_costbl[sc] = cos (PHI);
			sc++; if (phi<180) halfsino++;
		}
	}
	else
	{
		for (phi=0; phi<360; phi+=deltaphi)	
		sc++; if (phi<180) halfsino++;
	}


	line=0; sc=0;
#ifdef MIRRORME
	for (phi=0; phi<180; phi+=deltaphi)	/* Half revolution - use mirroring */
#else
	for (phi=0; phi<360; phi+=deltaphi)	/* Full revolution */
#endif
	{

		dp (1, "Now processing projection %d at phi=%f\n",sc+1,phi);
		for (j=0; j<arraywidth; j++)	/* One scan along the detectors */
		{

			if (!macrorun) indicate_progress (phi*j/(360.0*arraywidth));

			/* This is the formal way, sum Wij * fi,
				which, of course, is highly inefficient, but
				this is the way it's defined... */
			
			sum=0; cnt=0;
			for (y=0; y<img->ymax; y++)
				for (x=0; x<img->xmax; x++)
					if (in_recon_circle (x-img->xmax/2,y-img->ymax/2,arraywidth))
					{
						buf = Wij2 (x,y,j,sc,arraywidth);
						if (buf!=0.0)
						{
							sum += buf * readbuf_flt(img,x,y,z);
							cnt++;
						}
					}

			p[j + arraywidth*line] = sum;
#ifdef MIRRORME
			p[ (arraywidth-j) + arraywidth*(line+halfsino)] = sum;	/* Use symmetry */
#endif	
		}

		line++; sc++;
		if (cnt > arraywidth) dp (1, "Warning: Too many Wij at detector %d, angle %f\n",j,phi);
	}

	if (myjob)
	{
		free (art_sintbl); art_sintbl=NULL;
		free (art_costbl); art_costbl=NULL;
	}

	if (!macrorun) reset_progress ();

	handback (img, &sinogram);
#undef MIRRORME
}



/* Reconstruction by ART and SART. Upon call, img contains
	a sinogram (only the upper half is used, though). Upon exit,
	the reconstructed plane is passed back in img.
	Whichart determines the correction technique.
		Whichart==0 indicates SART, and whichart==1 indicates Kaczmarz ART.
		Whichart==2 is for SART with Hamming window for backprojection
	Halfview determines the type of sinogram:
		Halfview==1 means the sinogram has 180° coverage, use all lines.
		Halfview==0 means the sinogram has 360° coverage, use only the upper half.
	*/


void art_recon (image_type *img, int whichart, int halfview)
{
image_type recon;
int x,y,xi,yi,xmax,ymax,itercnt;
int i,j,sc,halfsine,cnt,c,skew;
AFLOAT anginc,buf,wij,sum,wijsum, phi, epsilon, nepsilon, lambda, span;
AFLOAT H,s,l,x0,y0;		/* For Hamming */
AFLOAT angcov;			/* Angular coverage, 360° or 180° */
AFLOAT *p, *haux;
float *q;
float xdata[30],ydata[30];
gchar sbuf[256];


	xmax = ymax = img->xmax;	/* source sinogram */
	if (halfview)
	{
		angcov = 180.0;
		anginc = angcov/img->ymax;
		skew = MAX (1, (int)(img->ymax / 2.44) );
	}
	else
	{
		angcov = 360.0;
		anginc = angcov/img->ymax;
		skew = MAX (1, (int)(img->ymax / 4.88) );
	}

	skew = find_prime (skew);	/* We need a prime number or it won't always cycle through all projections */
	lambda=0.7;					/* Undercorrection factor, 0<lambda <=1 */

	/* Allocate memory for the correction and reconstruction images */

	if (allocate_image (&recon, B_FLOAT, xmax,ymax,1)==-1)
		fatal_memerr();
	q = recon.data;
	
	/* initialize the recon image with average density */
	
	buf=0; cnt=0;
	for (y=0; y<img->ymax*img->xmax; y++)
	{
		if ( ((float*)img->data)[y] > -1e-8)
		{
			buf += ((float*)img->data)[y]; cnt++;
		}
	}
	if (cnt>0)
	{
		buf /= (SQR(cnt) / img->ymax);

		for (y=0; y<recon.ymax; y++)
			for (x=0; x<recon.xmax; x++)
			{
				if (in_recon_circle (x - recon.xmax/2, y - recon.ymax/2, recon.xmax))
					q[x + y*recon.xmax]=buf;
			}
	}


	if (art_use_refcorr==2) refraction_trick (img, 1);	/* Interpolation over shadow zones */
	
	p = calloc (xmax, sizeof (AFLOAT));		/* One line array worth of placeholders */
	if (!p) fatal_memerr();
	haux = calloc (3*xmax, sizeof (AFLOAT));	/* Auxiliary data for Hamming computation */
	if (!haux) fatal_memerr();

	/* Prepare the tables of pre-computed data */

	sc=0; halfsine=0;
	art_sintbl = calloc (img->ymax*4, sizeof(AFLOAT));
	art_costbl = calloc (img->ymax*4, sizeof(AFLOAT));
	for (phi=0; phi<360; phi+=anginc)	/* One full revolution including redundancy */
	{
		art_sintbl[sc] = SIN (phi);
		art_costbl[sc] = COS (phi);
		sc++; 
		if (phi<180) halfsine++;		/* Dummy op to determine the last line in sinogram */
	}

	j = precompute_art_area (0.5, 0.5);		/* Circle/beam overlap areas */

	if (halfview)
		halfsine = img->ymax;			/* Use the full sinogram with 180° coverage */
	else
		halfsine = img->ymax / 2;		/* Use only the upper half of the sinogram */


	itercnt=0;
	do				/* Start iterative reconstruction */
	{
	
		/* Note (1): A multithreaded version may project angles sc and sc+1 and then reconcile
			the projections by adding the results
			Note (2): SART is not truly SART as defined by Kak & Slaney. Here, only one projection
			is backprojected, whereas Kak & Slaney propose to compute the error for *all*
			projections (in a separate image) and THEN add the error to the reconstruction space.
			Note (3): Ideally, we should skew the reconstruction so that successive projections
			whose error is backprojected are at an angle of 73.8° (increment of sc is 73.8/anginc)
			Note (4): The underrelaxation lambda could be made adptive, starting with a low value
			and progressively increasing to 1.00 as the maximum error decreases
		*/

		sc=0; epsilon=0;
		for (c=0; c<halfsine; c++)		/* Run over all projection angles... */
		{
			/* Skew the projection number so that the angle increment is near 73.8 degrees.
				Conventional: use sc=c (angle increment of anginc)
			*/

			sc = (sc+skew) % halfsine;

			if (!macrorun) indicate_progress ((float)c/(float)halfsine);

			/* Compute the hamming window ray geometry */

			for (j=0; j<xmax; j++)	/* One scan along the detectors */
			{
				if (j >= xmax/2) 
					s = - (0.5+j-xmax/2);
				else
					s = 0.5+xmax/2-1-j;
				l = 2*sqrt(SQR(xmax/2)-SQR(s));		/* Length of ray inside recon window */
				x0 = s*art_costbl[sc];				/* Point on ray closest to origin */
				y0 = s*art_sintbl[sc];
				haux[3*j]=l; haux[3*j+1]=x0; haux[3*j+2]=y0;
			}

			dp (1, "Now processing projection %d at phi=%.2f\n",sc+1,sc*anginc);

			if ((whichart==0) || (whichart==2)) 		/* Perform SART or SART with Hamming */
			{

				/* Compute the forward projection difference for one ray into p[] */

				for (j=0; j<xmax; j++)	/* One scan along the detectors */
				{
					p[j] =  0.0;
					if ((readbuf_flt(img,j,sc,0) > -1e-8)	/* Don't use this projection if negative (=illegal) abs value */
					   || (art_use_refcorr==0))
					{
						sum=0; cnt=0; wijsum=0;
						for (y=0; y<recon.ymax; y++)
							for (x=0; x<recon.xmax; x++)
								if (in_recon_circle (x-recon.xmax/2,y-recon.ymax/2,recon.xmax))
								{
									wij = Wij_sched (x,y,j,sc,xmax);
									if (!F_ISZERO(wij))
									{
										sum += wij * q[x+recon.xmax*y];
										cnt++; wijsum += wij;
									}
								}
						if (wijsum>1e-8) 
							p[j] = (readbuf_flt(img,j,sc,0) - sum) / wijsum; 
					}
				}

				/* Backproject this correction over the reconstruction space */

				for (y=0; y<recon.ymax; y++)
					for (x=0; x<recon.xmax; x++)
						if (in_recon_circle (x-recon.xmax/2,y-recon.ymax/2,recon.xmax))
						{
							sum=0; wijsum=0; xi=x-xmax/2; yi=y-ymax/2;
							for (j=0; j<xmax; j++)
								if (!F_ISZERO (p[j]))	/* Do not go through this trouble if we backproject a zero */
								{
									wij = Wij_sched (x, y, j, sc, xmax);
									if (fabs(wij) > 1e-8)
									{
										if (whichart==2)		/* With Hamming window */
										{
											H = hamming (sqrt ( SQR( xi-haux[j*3+1] ) + SQR( yi-haux[j*3+2] ) ), haux[j*3]);
											wijsum += H*wij;
											sum += H*wij*p[j];
										}
										else					/* Without Hamming window */
										{
											wijsum += wij;
											sum += wij*p[j];
										}
									}
								}

							if (fabs(wijsum) >1e-8) 
							{
								sum = lambda*sum/wijsum;
								q[x + recon.xmax*y] += sum;
								epsilon = MAX(epsilon,sum);	/* Monitor convergence */
							}
						}
			}
			else if (whichart==1)			/* perform Kaczmarz ART correction */
			{
				for (j=0; j<xmax; j++)	/* One scan along the detectors */
				{
					if (readbuf_flt(img,j,sc,0) > -1e-8)	/* Don't use this projection if negative (=illegal) abs value */
					{
						sum=0; cnt=0; wijsum=0;
						for (y=0; y<recon.ymax; y++)
							for (x=0; x<recon.xmax; x++)
							if (in_recon_circle (x-recon.xmax/2,y-recon.ymax/2,recon.xmax))
								{
									wij = Wij_sched (x,y,j,sc,xmax);
									if (wij>1e-8)
									{
										sum += wij * q[x+recon.xmax*y];
										cnt++; wijsum += SQR(wij);
									}
								}

						if (fabs(wijsum)>1e-8)
						{ 
							buf = lambda*(readbuf_flt (img,j,sc,0) - sum)/wijsum;

							/* For each detector *AND* each angle, distribute over recon image */

							for (y=0; y<recon.ymax; y++)
								for (x=0; x<recon.xmax; x++)
									if (in_recon_circle (x-recon.xmax/2,y-recon.ymax/2,recon.xmax))
									{
										sum = buf*Wij_sched(x,y,j,sc,xmax);
										q[x+recon.xmax*y] += sum;
										epsilon = MAX(epsilon,sum);
									}
						}
					}
					
				}

			}
		}

		if (art_use_smtrick) art_tricks (&recon, 1);		/* Smoothing trick */
		if (art_use_ctrick ) art_tricks (&recon, 2);		/* Constraint trick */

		/* Normalized error and adaptive underrelaxation coefficient */

		make_minmax (&recon);
		span = recon.maxval.fmax - recon.minval.fmin;
		nepsilon = 100.0*epsilon/span;
		lambda = adapt_underrelaxation (lambda, nepsilon);

		if (!macrorun)
		{
			sprintf (sbuf, "Iteration %d, epsilon=%2.7f", itercnt, epsilon);
			statusbar_display (sbuf);
			secondary_display (&recon, 0, "Reconstruction progress");
			xdata[itercnt]=itercnt; ydata[itercnt]=epsilon;
		}

		if (!macrorun) reset_progress ();

		itercnt++;
		dp (1, "Finished iteration %d. Epsilon is %f (%5.1f percent of span). Underrelaxation %5.2f\n",itercnt, epsilon, nepsilon, lambda);



	}
	while ((itercnt<art_maxiter) && (epsilon > span/100.0));

	if (!macrorun) graphme (xdata,ydata,itercnt, "Iteration", "Epsilon", "Convergence of SART Reconstruction", PLOT_LOG_Y);	
	free (art_sintbl); art_sintbl=NULL;
	free (art_costbl); art_costbl=NULL;
	free (p); free (haux);
	handback (img, &recon);
}


void wijtest (image_type *img, int numdet, int numang, int whichart)
{
image_type recon;
int x,y,xmax,ymax,N;
int i,j,sc,halfsine,cnt;
AFLOAT anginc,wij,sum,wijsum, phi;
AFLOAT s,l, x0, y0;
float *q;
float xdata[30], ydata[30];




	xmax = ymax = numdet;	/* source sinogram */
	N = numdet/2;			/* Center detector number */
	anginc = 360.0/numang;

	/* Allocate memory for the correction and reconstruction images */

	if (allocate_image (&recon, B_FLOAT, xmax,ymax,1)==-1)
		fatal_memerr();
	q = recon.data;

	sc=0; halfsine=0;
	art_sintbl = calloc (img->ymax, sizeof(AFLOAT));
	art_costbl = calloc (img->ymax, sizeof(AFLOAT));
	for (phi=0; phi<360; phi+=anginc)	/* One full revolution including redundancy */
	{
		art_sintbl[sc] = sin (PHI);
		art_costbl[sc] = cos (PHI);
		sc++; if (phi<180) halfsine++;
	}

	j = precompute_art_area (0.5, 0.5);		/* Circle/beam overlap areas */
	for (i=0; i<=j; i++)
	{
		xdata[i] = i*0.05;
		ydata[i] = art_overlap[i];
	}
	graphme (xdata, ydata, j+1, "Radius fraction", "Area", "Circle-Beam Overlap Area", 0);


	for (sc=0; sc<halfsine; sc++) 
	{
		for (j=0; j<xmax; j++)	/* One scan along the detectors */
		{
			if (j >= N) 
				s = - (0.5+j-N);
			else
				s = 0.5+N-1-j;
			l = 2*sqrt(N*N-s*s);
			x0 = s*art_costbl[sc];
			y0 = s*art_sintbl[sc];

			dp (1, "Ray j=%3d. Dist s=%6.2f, length l=%6.2f, origin (%6.2f,%6.2f)\n",j,s,l,x0,y0);

			sum=0; cnt=0; wijsum=0;
			for (y=0; y<recon.ymax; y++)
				for (x=0; x<recon.xmax; x++)
					if (in_recon_circle (x-recon.xmax/2,y-recon.ymax/2,recon.xmax))
					{
						wij = Wij_sched (x,y,j,sc,xmax);
						/*H = hamming ( sqrt ( SQR( (x-N)-x0 ) + SQR( (y-N)-y0 ) ), l);*/
						q[x + recon.xmax*y] += wij;
					}
		}
	}

	free (art_sintbl); art_sintbl=NULL;
	free (art_costbl); art_costbl=NULL;
	handback (img, &recon);
	return;
}



/********************************************************************************************

 ####    #####    ##*     ##*    #   #    #*#    #####   ####    #   #    ##*    #####  
 #   #   #       #   #   #   #   #   #   #   #     #     #   #   #   #   #   #     #    
 #   #   #       #       #   #   ##  #   #         #     #   #   #   #   #         #    
 ####    ####    #       #   #   # # #    #*#      #     ####    #   #   #         #    
 # #     #       #       #   #   #  ##       #     #     # #     #   #   #         #    
 #  #    #       #   #   #   #   #   #   #   #     #     #  #    #   #   #   #     #    
 #   #   #####    ##*     ##*    #   #    #*#      #     #   #    #*#     ##*      #    


	FILTERED BACKPROJECTION and FOURIER RECONSTRUCTION

********************************************************************************************/


/*	LIST OF RECONSTRUCTION FUNCTIONS

	shepplogan				Return a value of the Shepp-Logan deconvolution kernel
	ramachandran_lakshminarayanan		Ditto, but Ram-Lak
	reconstruct_prefilter	Filter the sinogram with one of the above
	ct_ammse_1d				1D AMMSE sinogram denoising (experimental)
	ct_bilateral_1d			Nonadaptive bilateral filter for sinogram denoising (experimental)
	ct_diffusion_1d			anisotropic diffusion for sinogram denoising (experimental)
	ctprepare				Convert intensity to absorbance, includes noise filters and centering
	reconstruct_slice		One-slice Filtered Backprojection algorithm
	fft_backproject			One-slice reconstruction via Fourier slice theorem
	fast_fft_preview		FFT preview for real-time previewing
*/





double shepplogan (int x)		/* Shepp-Logan deconvolution kernel */
{
	return -2.0/(SQR(PI) * (4.0*x*x-1.0));
}

double ramachandran_lakshminarayanan (int x)	/* Deconvolution kernel by Ramachandran and Lakshminarayanan */
{

	if (x==0)
		return 0.25;
	else if (x & 0x01)		/* odd x */
		return -1.0/SQR(PI*x);
	else
		return 0.0;			/* even x */
}


void reconstruct_prefilter (image_type* img, int xmax, int filter)
{
int x,y,z,q,i;
double sum;
image_type buf;
float* p;


	if (filter==0) return;		/* the no-op filter */

	if (allocate_image (&buf,B_FLOAT,img->xmax,img->ymax,img->zmax)==-1)
		fatal_memerr();
	p = buf.data;

	if (img->xmax & 0x01)		/* odd number of detectors */
		q = (img->xmax-1)/2;
	else						/* even number of detectors */
		q = img->xmax/2;
		
	/* It's a 1D convolution of the detector array line, line by line */

	for (z=0; z<img->zmax; z++)				/* Run over all slices */
	{
		for (y=0; y<img->ymax; y++)
		{
			for (x=0; x<img->xmax; x++)	/* convolve this line */
			{
				sum = 0;
				
				if (filter==1)		/* Shepp-Logan */		
					for (i=0; i<img->xmax; i++)
						sum += shepplogan (i-x)*readbuf_flt(img,i,y,z);
				else if (filter==2)	/* Ramachandran */
					for (i=0; i<img->xmax; i++)
						sum += ramachandran_lakshminarayanan (i-x)*readbuf_flt(img,i,y,z);
				
				p[x+img->xmax*(y + img->ymax*z)] = sum;
			}
		}
	}

	/* Hand it back */
	
	handback (img, &buf);
}



/* Sinogram denoising through adaptive minimal mean-squared-error filtering (AMMSE).
	This filter is a 1D version of the 2D version in the adaptive filters
	TODO: Th cimage 1.0.0, the AMMSE was massively overhauled. Those modifications must
	be reflected here, too, or this AMMSE will be widely worthless. */

void ct_ammse_1d (float *p, int N)
{
float gauss2[6] = {1.0, 0.882, 0.607, 0.325, 0.135, 0.044};		/* Gaussian, sigma = 2, l=11 */
float Sx, Sxx, w, minvar, maxvar, wsum, buf;
int cnt, i, j, k, ks;
float *locvar, *smoothed;			/* Since these convolutions are needed anyway, store & reuse them */


	locvar = calloc (N, sizeof (float));
	smoothed = calloc (N, sizeof (float));

	w=0; ks=5;						/* ks = kernel support */
	for (i=-ks; i<=ks; i++) w+=gauss2[abs(i)];

	/* In one pass, compute local variance, minimum local variance, and the smoothed version */

	for (j=0; j<N; j++)
	{

		Sx=0.0; Sxx=0.0; wsum=0.0; cnt=0;
		for (i=-ks; i<=ks; i++)
		{
			k = j+i; if (k<0) k=-k; else if (k>=N) k=2*N-k-1;		/* Mirror-boundary condition */
			Sx += p[k]; Sxx += SQR(p[k]); wsum += p[k]*gauss2[abs(i)]; cnt++;
		}

		smoothed[j] = wsum/w;
		locvar[j] = SDEV(Sx, Sxx, cnt);

		if (j==0)
		{
			maxvar = minvar = locvar[j];
		}
		else
		{
			if (minvar > locvar[j]) minvar=locvar[j];
			if (maxvar < locvar[j]) maxvar=locvar[j];
		}
	}

	dp (3, "Line: Min variance %f, max variance %f\n", minvar, maxvar);

	/* Pass 2: Return the AMMSE filtered data in place */

	if (!F_ISZERO(minvar))			/* No variance? Then this projection is virtually noise-free */
	{
		for (j=0; j<N; j++)
		{
			buf = p[j] - SQR(minvar/locvar[j]) * (p[j] - smoothed[j]);
			p[j] = buf;
		}
	}

	free (locvar);
	free (smoothed);

}


/* Sinogram denoising through bilateral filtering.
	This filter is a nonadaptive 1D version of the 2D version 
	in the adaptiver filters */

void ct_bilateral_1d (float *p, int N, float sig_d, float sig_r)
{
float w, wsum, buf;
float var_r, var_d;								/* domain and range filter variance times 2 */
int cnt, i, j, k, ks;
float *smoothed;


	smoothed = calloc (N, sizeof (float));
	ks = 5;										/* Kernel support, see AMMSE above. */
	var_d = 2.0*SQR(sig_d);
	var_r = 2.0*SQR(sig_r);

	/* The bilateral filter is a one-pass convolution filter with a data-dependent Gaussian */

	for (j=0; j<N; j++)
	{

		w=0.0; wsum=0.0; cnt=0;
		for (i=-ks; i<=ks; i++)
		{
			k = j+i; if (k<0) k=-k; else if (k>=N) k=2*N-k-1;		/* Mirror-boundary condition */
			buf = exp (-SQR(i)/var_d - SQR(p[j]-p[k])/var_r);		/* Kernel component */
			w += buf; wsum += p[k]*buf; cnt++;
		}

		smoothed[j] = wsum/w;
	}

	for (j=0; j<N; j++)
		p[j] = smoothed[j];

	free (smoothed);
}



/* Sinogram denoising through anisotropic diffusion
	This filter is a 1D version of the 2D version 
	in the adaptiver filters */

float diffusion_g (float gradi, float k)
{
float f;

	f = 1.0/ (1 + SQR(gradi/k));
	return f;
}



void ct_diffusion_1d (float *p, int N, float K, int iter)
{
float ir, il, cr, cl, buf, deltat;
int i,k;
float *smoothed;


	smoothed = calloc (N, sizeof (float));
	deltat = 0.2;

	for (k=0; k<iter; k++)
	{
		for (i=0; i<N; i++)
		{

			/* Compute the two intensity gradients (note sign difference to paper */

			buf = p[i];
			ir = il = 0.0; 
			if (i>0) il = buf - p[i-1];
			if (i<N-1) ir = buf - p[i+1];

			/* Compute the four conduction factors - those reaching out of the image are zero */

			cr = cl = 0.0;
			if (i>0) cl = anisotropic_diffusion_g (il,K);
			if (i<N-1) cr = anisotropic_diffusion_g (ir,K);

			/* Compute the new iteration value at (x,y) and store in target.
					Note sign difference to paper */

			smoothed[i] = buf - deltat*(cl*il + cr*ir);
		}

		for (i=0; i<N; i++)
			p[i] = smoothed[i];

	}

	free (smoothed);
}




/* CT-Prepare works on the individual projections of a sinogram.
   Preparation of real-world data is performed as follows:

	1. If requested, perform sinogram-level denoising (first step, because detector noise is additive)
	2. If this is not data from a log-amplifier, take the logarithm of all data
	3. Both margins are averaged to get the value of air. Then each pixel
   		is subtracted from the air value (under the assumption that we have 
   		log data this is equivalent to log (I0/I)).
	4. Center the sinogram. Compute the center of gravity of all views
	   and shift it into the center
  
	- OR -

	5. Determine the center of gravity for each line and fit it to a sine wave

   We do not assume img to be of B_FLOAT type but upon output it is.
*/


void ctprepare (image_type* img, int is_logdata, int shift, int manualshift, float pixsize, int sino_denoise)
{
int i,x,y,z,cl,cr,cnt, delta;
image_type I;
float l,r,buf, air, min,max, cogx,cogy, thresh, pixfac;
float *p, *p1;
long idx;


	if (allocate_image (&I,B_FLOAT, img->xmax, img->ymax, img->zmax)==-1)
		fatal_memerr();
	p1 = calloc (img->xmax, sizeof (float));
	if (p1==0) fatal_memerr();

	pixfac = 10.0/pixsize;	/* So we get mu in inverse centimeters */


	/* In preparation for filtering and conversion to absorbance,
		normalize the data. We eventually need to have a projection P
		as P(t) = -ln[I(t)/I0] = ln[I0] - ln[I(t)],
		so we begin by creating the air-normalized intensity I(t)/I0,
		which ranges from 0 (photon-starved) to 1 (air).
		We directly use the new image placeholder I, because it is float
		anyway.
	*/

	p = I.data;						/* Result data */
	for (z=0; z<img->zmax; z++)		/* Slice by slice */
	{
		for (y=0; y<img->ymax; y++)		/* View by view */
		{
			l=0.0; r=0.0; cl=0; cr=0;
			idx = img->xmax*(y+img->ymax*z);
			for (x=0; x<img->xmax; x++)		/* average 10 pixels each side */
			{
				buf = readbuf_flt(img,x,y,z);		/* Optimize...? */
				if (x<10) { l+=buf; cl++; }
				if (x>img->xmax-11) { r+=buf; cr++; }
			}

			/* We could linearly interpolate between l and r, but...  */

			air = 0.5*(l/cl+r/cr);
			for (x=0; x<img->xmax; x++)
				p[x+idx] = readbuf_flt(img,x,y,z)/air;
		}
	}


	/* If requested, perform denoising BEFORE any other operation.
		We can do the filtering in-place in I */

	p = I.data;
	if (sino_denoise)
	{
		for (z=0; z<img->zmax; z++)		/* Slice by slice */
		{
			for (y=0; y<img->ymax; y++)		/* and line by line */
			{
				if (sino_denoise==1)
					ct_ammse_1d (p, img->xmax);
				else if (sino_denoise==2)
					ct_bilateral_1d (p, img->xmax, 3.0 /* domain sigma */, 2.0 /* range sigma */);
				else if (sino_denoise==3)
					wvlt_1D_denoising (p, img->xmax, 2 /* Daub-8 */, 3 /* BayesShrink hard */, 1.0);
				p+=img->xmax;
			}
		}
	}


	p = I.data;						/* Rewind the pointer as we use readbuf to access data in I */

	for (z=0; z<img->zmax; z++)		/* Slice by slice */
	{

		/* CONVERT TO ABSORBANCE. Note that we'd theoretically have to
			clamp values where I(t)/I0 > 1, because this would indicate negative
			absorbance. Well, let's do it anyhow. */

		for (y=0; y<img->ymax; y++)			/* View by view */
		{
			idx = img->xmax*(y+img->ymax*z);
			for (x=0; x<img->xmax; x++)		/* and pixel by pixel */
			{
				if (is_logdata)
					buf = (1.0 - p[x+idx])*pixfac;
				else
				{
					p[x+idx] = -pixfac*SAFELOG(p[x+idx]);
					if (p[x+idx] < 0) p[x+idx]=0;				/* clamp */
				}
			}
		
		}

		/* We need a separate pass to obtain min/max because we need p[0] */

		min=p[0]; max=min;
		for (y=0; y<img->ymax; y++)		/* View by view */
		{
			idx = img->xmax*(y+img->ymax*z);
			for (x=0; x<img->xmax; x++)
			{
				buf = p[x+idx];
				if (min>buf) min=buf;
				if (max<buf) max=buf;
			}
		}

		/* Shift the COG into the center */

		if (manualshift)
		{
			delta = shift;
			dp (1, "Sinogram: Manual shift by %d pixels\n",delta);
		}
		else
		{
			cogx=0; cogy=0; cnt=0;
			thresh = 0.25;				/* Threshold: Percent of the span to be considered object */
			for (y=0; y<img->ymax; y++)
			{
				idx = img->xmax*(y+img->ymax*z);

				x=0; i=0;
				while ((i<3) && (x<img->xmax))		/* approach from the left */
				{
					buf = p[x+idx];
					if (buf > thresh*(max-min)) 
						i++;
					else
						p[x+idx]=0;	/* threshold */
					x++;
				}
				l = x-2;
				x=img->xmax-1; i=0;
				while ((i<3) && (x>0))		/* approach from the right */
				{
					buf = p[x+idx];
					if (buf > thresh*(max-min)) 
						i++;
					else
						p[x+idx]=0;	/* threshold */
					x--;
				}
				r = x+2;

				for (x=l; x<=r; x++) { cogx+=x; cogy+=y; cnt++; }

			}

			cogx /= cnt; cogy/=cnt;
			delta = (int)(0.5+(img->xmax/2.0) - cogx);
			if ((delta > 256) || (delta < -256)) delta=0;	/* ignore obvious trash values */
			dp (1, "Sinogram: found centroid at %6.1f,%6.1f. Shift by %d pixels\n",cogx,cogy,delta);
		}

		/* Shift sinogram by delta */

		if (delta!=0)
		{
			for (y=0; y<img->ymax; y++)
			{
				idx = img->xmax*(y+img->ymax*z);
				for (x=0; x<img->xmax; x++) p1[x]=0.0;
				for (x=0; x<img->xmax; x++)
				{
					if ((x+delta>=0) && (x+delta <img->xmax))
						p1[x+delta] = p[x+idx];
				}
				for (x=0; x<img->xmax; x++) p[x+idx]=p1[x];
			}
		}

	}
	handback (img, &I);

}


/* The filtered backprojection algorithm */

/*	We are using the backprojection approach in the
	spatial domain. There are two steps required,
	the pre-filtering with a deconvolution kernel, e.g.
	the Shepp-Logan kernel; and the actual backprojection,
	which maps the sinogram s(d,theta) to f(x,y).
	The deconvolution prefilter accepts a sinogram and
	returns the filtered sinogram in the same pointer.
	The backprojection algorithm accepts a sinogram
	and returns the reconstructed slice.
*/


void reconstruct_slice (image_type* sino, int xmax, int halfview)
{
int x,y,z,i,zmax;
long idx;
char err;
float* p;	/* the reconstructed image */
image_type recon;
double sum, deltaphi,phi,t,a,fx,fy, xmaxhalf;
double *sintbl, *costbl;

	/* Get memory for reconstructed image, prepare some constants */

	zmax = sino->zmax;
	if (allocate_image (&recon, B_FLOAT, xmax, xmax, zmax)==-1)	/* square image */
		fatal_memerr();
	p = recon.data;
	if (halfview)
		deltaphi = PI/sino->ymax;	/* Angle increment between views for 180 degree coverage */
	else
		deltaphi = 2*PI/sino->ymax;	/* Angle increment between views for 360 degree coverage */
		
	xmaxhalf = 0.5*sino->xmax;
	a = xmaxhalf/recon.xmax; /* Expansion factor, constant within the loops */
	a = sino->xmax/recon.xmax; /* Expansion factor, constant within the loops */

	/* For faster computation, we prepare sin/cos tables in advance.
		In the inner loop, phi is i*deltaphi, so each table element
		at [i] contains the value sin(phi*deltaphi) */
	
	sintbl = calloc (sino->ymax, sizeof(double));
	costbl = calloc (sino->ymax, sizeof(double));
	for (i=0; i<sino->ymax; i++)
	{
		sintbl[i] = sin (i*deltaphi);
		costbl[i] = cos (i*deltaphi);
	}


	/* the actual backprojection sums sinogram detectors over all
	   views. The detector used depends on the coordinate */

	for (z=0; z<zmax; z++)
	{
		for (y=0; y<recon.ymax; y++)
		{
			if (!macrorun) indicate_progress ( ((float)y + (float)z*recon.ymax)/((float)recon.ymax*(float)zmax) );

			for (x=0; x<recon.xmax; x++)
			{
				sum = 0; phi=0; err=0;
				fx = (float)x - 0.5*recon.xmax;	/* center within sinogram */
				fy = (float)y - 0.5*recon.ymax;
				for (i=0; i<sino->ymax; i++)	/* sum over all views */
				{
					phi += deltaphi; /* phi = i*deltaphi; use tables instead */
					t = (fx*costbl[i] + fy*sintbl[i])*a; /* detector */
					if ((t<-xmaxhalf) || (t>=xmaxhalf)) 
					{
						dp (3, "Warning: Illegal detector %5.1f at (%d,%d), view %d\n",t,x,y,i);
						err=1;
					}
					else
						sum += ireadbuf(sino,t+xmaxhalf,i,z,intp_pref);
				}
				/* sum *= 2*PI/sino->ymax; according to Natter - The 2 PI seems to be wrong*/

				/* Somehow, X and Y seem to be flipped... Quick fix here */

				idx = (recon.xmax-x-1) + recon.xmax*(recon.ymax-y-1) + recon.xmax*recon.ymax*z;
				if (!err)
					p[idx] = PI*sum/sino->ymax;
				else
					p[idx] = 0;
			}

		}
	}

	if (!macrorun) reset_progress ();

	/* Hand it back */
	
	free (sintbl);
	free (costbl);
	handback (sino, &recon);

}



/* The alternative approach is the Radon-backprojection in
   Fourier space. The sinogram is transformed line-by-line.
   Then converted into a polar array and back-FFT'd.
 
*/


void fft_backproject (image_type* sino, int xmax, int halfview, int preview)
{
fftw_complex *p;		/* typed pointer to 2D fft data area */
fftw_complex *op;		/* typed pointer to 2D line-by-line FFT'd sinogram */
fftw_complex *lp;		/* FFTW input line */
fftw_complex *lpz;		/* typed pointer to fft line area */

fftw_plan plan1d;
fftw_plan plan2d;

int i, xm, xh,yh,x,y,x1,y1, np, th;
int z, zidx, zmax, yidx;
float f,a,r;
image_type recon; float* fp;
image_type dbgimg; float* dbgp;


	xm = sino->xmax;			/* Number of detectors */
	zmax = sino->zmax;			/* Number of slices */
	if (preview) zmax=1;		/* Only slice 0 allowed in preview mode */
	if (halfview)
		np = sino->ymax;	/* Number of projections for 180 degree coverage, uses all projections */
	else
		np = sino->ymax/2;	/* Number of projections for 360 degree coverage, uses upper half */
		
	if (xm & 1) xm-=1;	/* if odd number, simply ignore one detector */
	xh = xm/2; yh = xh;

	/* Get the memory we need */
	
	p = calloc (SQR(xm), sizeof (fftw_complex));	/* For the 2D transformation */
	if (p==0) fatal_memerr();

	lp = calloc (xm*2, sizeof (fftw_complex));		/* For the 1D transform */
	lpz= calloc (xm*2, sizeof (fftw_complex));		/* For the 1D transformed */
	op = calloc (xm*np*2, sizeof (fftw_complex));		/* To store the 1D transform */
	if ((lp==0) || (lpz==0) || (op==0)) fatal_memerr();

	fp = NULL; recon.data = NULL;
	if (allocate_image (&dbgimg,B_FLOAT,xm,xm,2)==-1)
		fatal_memerr();
	dbgp = dbgimg.data;

	if (allocate_image (&recon,B_FLOAT,xm,xm,zmax)==-1)
		fatal_memerr();

	/* Estimate transformation strategies.
	   For the one-dimensional forward transformation we use 1D 
	   real-to-complex from lp into lpz.
	   For the two-dimensional inverse transform we use 2D
	   in-place transform within p */
	
	
	plan1d = fftw_plan_dft_1d (xm, lp, lpz, FFTW_FORWARD, FFTW_ESTIMATE);

	plan2d = fftw_plan_dft_2d (xm,xm, p, p, FFTW_BACKWARD, FFTW_ESTIMATE);

	/* Begin reconstruction slice-by-slice */

	for (z=0; z<zmax; z++)
	{
		dp (1, "Radon: Reconstructing slice %d\n", z);


		/* Fill the line array. This is needed because "sino" may have arbitrary type.
		Note that the original sinogram center needs to be moved to the margins -
		I have no idea why, but this is how it works. Also, we don't need a full
		360 degree sinogram, so crop the point-mirrored bottom part
		(this is why np = ymax/2 when halfview==0) */

		dp (2, "Radon: Filling the line array ");
		for (y=0; y<np; y++)
		{
			for (x=0; x<xm; x++)
			{
			
				if (x<xh)
				{
					c_re(lp[x+xh]) = readbuf_flt (sino,x,y,z);		/* Real component only */
					c_im(lp[x+xh]) = 0.0;
				}
				else
				{
					c_re(lp[x-xh]) = readbuf_flt (sino,x,y,z);
					c_im(lp[x-xh]) = 0.0;
				}
			}

			fftw_execute (plan1d);
			
			/* Now shuffle the transform so that the zero frequency is in the center.
			The actual zero element is therefore at x=xh in each line.
			In addition, multiply with a generalized Hamming window
			(see G. Herman, 1980) */
			
			for (x=0; x<xm; x++)
			{
				/*
				lpz[x].re *= (0.8+0.2*cos(PI*(x-xh)/xm));
				lpz[x].im *= (0.8+0.2*cos(PI*(x-xh)/xm));
				*/
				if (x<xh)
				{
					c_re(op[x+xh+y*xm]) = c_re(lpz[x]);
					c_im(op[x+xh+y*xm]) = c_im(lpz[x]);
				}
				else
				{
					c_re(op[x-xh+y*xm]) = c_re(lpz[x]);
					c_im(op[x-xh+y*xm]) = c_im(lpz[x]);
				}
			}

		}
		dp (2,"\n");


		/* Now prepare the 2D array. Each coordinate (x,y) of the new array
		has a cooresponding pair (r,P) in the transformed sinogram, where
		x = r*cos P and y = r * sin P. P is the angle (corresponds to
		projection line) and r is the frequency index. The reason is that
		each transformed line of the sinogram corresponds to one line
		in the 2D array, at an angle. For example, the sinogram line
		at 45 deg projection angle is a line in the 2D array at 45 degrees.
		This means that each pixel in the 2D array has polar coordinates r,P
		where r is the frequency index and P the sinogram line */
		
		dp (2, "Radon: Preparing polar array ");
		xh = xm/2; yh=xh;

		/* Reverse approach: For each pixel in the 2D matrix
		find the closest match in the 1D transform. This is the preferable
		approach as it gets dramatically better results (empirical) */

		for (y=0; y<xm; y++)
			for (x=0; x<xm; x++)
			{
				y1=y-yh; x1=x-xh; 				/* Wrap: Center in the middle */
				if ((x1==0) && (y1==0))
				{
					r=0.0; th=0;
				}
				else
				{
					r = sqrt(x1*x1+y1*y1); 			     /* Make polar coords */
					th= (int) (0.5 + atan2(y1,x1)*np/PI);
				}

				/* atan2 returns a value between -PI and PI, so our view index th
					may be negative. In this case, r is in reality negative.
					Translate r so that r=0 indexes xh in the 1D transform */

				if (th<0)
				{
					th += np; r = xh-r;
				}
				else
					r = xh+r;

				if ((th>=0) && (th<np)  && (r>=0) && (r<xm))
				{
					/* We interpolate between detectors, but not between views */
					
					i = (int)floor(r); a=r-i;
					c_re(p[x + xm*y]) = (1-a)*c_re(op[i+xm*th]) + a*c_re(op[i+1+xm*th]);
					c_im(p[x + xm*y]) = (1-a)*c_im(op[i+xm*th]) + a*c_im(op[i+1+xm*th]);
				}
				else
				{
					c_re(p[x + xm*y]) = 0.0;
					c_im(p[x + xm*y]) = 0.0;
				}
			}
			
		/* To reduce numerical errors, copy the zero degree row verbatim */
		
		for (x=0; x<xm; x++)
		{
			c_re(p[x + xm*xh]) = c_re(op[x]);
			c_im(p[x + xm*xh]) = c_im(op[x]);
		}			

		dp (2, "OK\n");

		if (!(preview) && (!macrorun))
		{
			/* Before doing ifft, visualize the polar array. This should have the
			zero-frequency point in the center */

			for (y=0; y<xm; y++)
				for (x=0; x<xm; x++)
				{
					dbgp[x+xm*y] = c_re(p[x+xm*y]);
					dbgp[x+xm*y+xm*xm] = c_im(p[x+xm*y]);
				}
			make_minmax (&dbgimg);
			secondary_display (&dbgimg, 0, "2D Radon-transform");
			freebuf (&dbgimg);
		}


		/* The inverse FFT requires that we swap quadrants so that the
		zero-frequency point ends up in the (0,0) corner */

		dp (2, "...inverse fft...\n");
		
		for (y=0; y<yh; y++)
			for (x=0; x<xh; x++)
			{
				fftw_real R;
				
				R = c_re(p[x+xm*y]); c_re(p[x+xm*y]) = c_re(p[x+xh+xm*(y+yh)]); c_re(p[x+xh+xm*(y+yh)]) = R;
				R = c_im(p[x+xm*y]); c_im(p[x+xm*y]) = c_im(p[x+xh+xm*(y+yh)]); c_im(p[x+xh+xm*(y+yh)]) = R;
				
				R = c_re(p[x+xm*(y+yh)]); c_re(p[x+xm*(y+yh)]) = c_re(p[x+xh+xm*y]); c_re(p[x+xh+xm*y]) = R;
				R = c_im(p[x+xm*(y+yh)]); c_im(p[x+xm*(y+yh)]) = c_im(p[x+xh+xm*y]); c_im(p[x+xh+xm*y]) = R;
			}


		fftw_execute (plan2d);		/* Inverse 2D transform */
		dp (2, "Done\n");


		/* Finally, hand back via "recon" from the array at pointer p. Normalize the data.
			Note - for some reasons the image is rotated by 180 degs,
			so we do a quick fix here */
		
		/* Also note that this operation returns two slices, the real and the imaginary
			part. Ideally, the imaginary part would be zero.
			The commented-out sections would return a single slice with the magnitude. */
		
		fp = recon.data;		/* Target image for reconstructed slice */
		zidx = recon.xmax*recon.ymax*z;
		
		f = 1.0/(xm*xm);
		for (y=0; y<xm; y++)
		{
			if (y<yh) y1=yh+y; else y1=y-yh;
			yidx = xm*(recon.ymax-y1-1);

			for (x=0; x<xm; x++)
			{
				if (x<xh) x1=xh+x; else x1=x-xh;

				fp[(recon.xmax-x1-1)+yidx+zidx] = c_re(p[x+xm*y])*f;
				/*
				fp[(recon.xmax-x1-1)+yidx+xm*xm] = c_im(p[x+xm*y])*f;

				fp[(recon.xmax-x1-1)+xm*(recon.ymax-y1-1)] = sqrt(SQR(c_re(p[x+xm*y]))
																+SQR(c_im(p[x+xm*y])))*f;
				fp[(recon.xmax-x1-1)+xm*(recon.ymax-y1-1)+xm*xm] = 0;
				*/

			}
		}

	}				/* End for all slices */

	free (lp); free (lpz);		/* The 1D line arrays are no longer needed */
	fftw_destroy_plan (plan1d);
	fftw_destroy_plan (plan2d);
	free (op);


	handback (sino, &recon);
	make_minmax (sino);
	free (p);		
}


/* Fast preview using Radon backtransformation */


void fast_fft_preview (image_type* sino, int shift, int realworld, int islogdata, int halfview)
{
int x,y,e;

long idx1, idx2;
unsigned char *p1, *p2;
image_type smallsino;

	if (allocate_image (&smallsino, sino->imgtype, sino->xmax, sino->ymax,1)==-1)
		fatal_memerr();
	e = DATALEN(sino->imgtype);
	p1 = sino->data; p2 = smallsino.data;
	
	for (y=0; y<sino->ymax; y++)
		for (x=0; x<sino->xmax; x++)
		{
			idx1 = e*(x + sino->xmax*y);
			idx2 = x+shift; 
			if (idx2<0) idx2+=sino->xmax; else if (idx2>=sino->xmax) idx2-=sino->xmax;
			idx2=e*(idx2+sino->xmax*y);
			memcpy (&p2[idx2], &p1[idx1], e);
		}
	
	if (realworld) 
		ctprepare (&smallsino,islogdata,0,1,0.1,0);		/* Note that we ignore pixel scaling and denoising */

	fft_backproject (&smallsino, smallsino.xmax, halfview, 1);
	secondary_display (&smallsino, 0, "Reconstruction preview");
	freebuf (&smallsino);

}



/****************************************************

	Sinogram / CT reconstruction GUI
	
*****************************************************/


GtkAdjustment *sino_anginc_adj, *sino_numdets_adj, *recon_resolution_adj;
GtkAdjustment *recon_shift_adj, *recon_pixsize_adj;
GtkWidget *advanced_art_wnd, *art_opt_button;
int sino_alg, sino_fanbeam, recon_kernel, sino_realworld, sino_logdata, recon_manualshift, sino_straighten;
int sino_advart_exists, sino_halfview, sino_denoise;


void on_sino_cancelbutton_clicked  (GtkButton *button, gpointer user_data)
{

	gtk_widget_destroy (GTK_WIDGET (user_data));
}


void on_sino_okbutton_clicked  (GtkButton *button, gpointer user_data)
{
int numdets, anginc;

	undo_nexus (&mainimg);

	numdets = sino_numdets_adj->value;
	anginc = sino_anginc_adj->value;

	if (sino_alg==0)
		make_sinogram_simple (&mainimg, cur_slice, numdets, anginc, sino_fanbeam);
	else if (sino_alg==1)
		art_project (&mainimg, cur_slice, numdets, anginc);
	else if (sino_alg==2)
		make_sinogram_radon (&mainimg, cur_slice, numdets, anginc);
	else if (sino_alg==3)
		make_pet_sinogram (&mainimg, cur_slice, numdets, anginc);

	make_imlib_buffer_from_main_image (cur_slice);
	gtk_widget_destroy (GTK_WIDGET (user_data));
}


void on_sino_fanbeam_clicked  (GtkButton *button, gpointer user_data)
{
	sino_fanbeam = !sino_fanbeam;
}


void on_sino_realwld_clicked  (GtkButton *button, gpointer user_data)
{
	sino_realworld = !sino_realworld;
}


void on_sino_logdata_clicked  (GtkButton *button, gpointer user_data)
{
	sino_logdata = !sino_logdata;
}


void on_sino_straighten_clicked  (GtkButton *button, gpointer user_data)
{
	sino_straighten = !sino_straighten;
}


void set_sino_menu (GtkWidget *widget, gpointer user_data)
{

	sino_alg = *(int*) user_data;
	
}

void on_sino_halfview_clicked  (GtkButton *button, gpointer user_data)
{
	sino_halfview = !gtk_toggle_button_get_active ((GtkToggleButton*)user_data);
	
	if (sino_halfview)
		dp (2, "Sinogram is assumed to cover 180 degrees\n");
	else
		dp (2, "Sinogram is assumed to cover 360 degrees\n");
}



void create_sinogram_wnd ()
{
GtkWidget* wnd;
GtkWidget* vbox;
GtkWidget* table;
GtkWidget *label1, *label2, *label3, *label4;
GtkWidget *sino_anginc, *sino_numdets;
GtkWidget* fanbeam_checkbutton;
GtkWidget *sino_okbutton, *sino_cancelbutton;
GtkWidget *glade_menuitem;
GtkWidget *sino_menu, *sino_menu_menu;


	wnd = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	gtk_object_set_data (GTK_OBJECT (wnd), "wnd001", wnd);
	gtk_container_set_border_width (GTK_CONTAINER (wnd), 3);
	gtk_window_set_title (GTK_WINDOW (wnd), "Create Sinogram");
	gtk_window_set_modal (GTK_WINDOW (wnd), TRUE);
	vbox = vbox_in_window (wnd,"vbox001");
	table = MakeTable (wnd,vbox,"table001",2,5);
	label1 = MakeLabel (wnd,table,"Angular increment","label1",0,0);
	label2 = MakeLabel (wnd,table,"Number of detectors","label2",0,1);
	label3 = MakeLabel (wnd,table,"Algorithm","label3",0,2);
	label4 = MakeLabel (wnd,table,"Fan-beam geometry","label4",0,3);

	sino_anginc_adj = GTK_ADJUSTMENT (gtk_adjustment_new (5, 1, 180, 1, 10, 0));
	sino_anginc = MakeHScale (sino_anginc_adj,wnd,table,"sino_anginc",1,0);
	gtk_scale_set_digits (GTK_SCALE (sino_anginc), 0);
	
	sino_numdets_adj = GTK_ADJUSTMENT (gtk_adjustment_new (mainimg.xmax,1,1024, 1, 10, 0));
	sino_numdets = MakeHScale (sino_numdets_adj,wnd,table,"sino_numdets",1,1);
	gtk_scale_set_digits (GTK_SCALE (sino_numdets), 0);

	sino_menu = Make_Optionmenu (wnd, table,"sino_kernel_menu",1,2);
	gtk_container_set_border_width (GTK_CONTAINER (sino_menu), 1);
	sino_menu_menu = gtk_menu_new ();
	glade_menuitem = Option_menuitem (sino_menu_menu, GTK_SIGNAL_FUNC(set_sino_menu), "Raytrace (Slow)", (gpointer) &menuopt0);
	glade_menuitem = Option_menuitem (sino_menu_menu, GTK_SIGNAL_FUNC(set_sino_menu), "Arithmetic (Very slow)", (gpointer) &menuopt1);
	glade_menuitem = Option_menuitem (sino_menu_menu, GTK_SIGNAL_FUNC(set_sino_menu), "Radon (Fast)", (gpointer) &menuopt2);
	gtk_option_menu_set_menu (GTK_OPTION_MENU (sino_menu), sino_menu_menu);

	glade_menuitem = Option_menuitem (sino_menu_menu, GTK_SIGNAL_FUNC(set_sino_menu), "Random Emission", (gpointer) &menuopt3);
	gtk_option_menu_set_menu (GTK_OPTION_MENU (sino_menu), sino_menu_menu);
	fanbeam_checkbutton = Checkbutton_in_table (wnd,table," ","sino_chkbutton",1,3,on_sino_fanbeam_clicked);
	sino_okbutton = Button_in_table (wnd,table,"OK","sino_okbutton",0,4,NULL);
	sino_cancelbutton = Button_in_table (wnd,table,"Cancel","sino_cancelbutton",1,4,NULL);

	gtk_signal_connect (GTK_OBJECT (wnd), "destroy",
                      GTK_SIGNAL_FUNC (on_sino_cancelbutton_clicked), wnd);
	gtk_signal_connect (GTK_OBJECT (sino_cancelbutton), "clicked",
                      GTK_SIGNAL_FUNC (on_sino_cancelbutton_clicked), wnd);
	gtk_signal_connect (GTK_OBJECT (sino_okbutton), "clicked",
                      GTK_SIGNAL_FUNC (on_sino_okbutton_clicked), wnd);

	gtk_widget_show (wnd);
}


/*************************************************************************************/




/* Advanced ART options (separate window) */


void destroy_advart_wnd (GtkWidget *widget, gpointer user_data)
{
	sino_advart_exists = -1;
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (art_opt_button), 0);
	sino_advart_exists = 0;
	gtk_widget_destroy ( GTK_WIDGET (advanced_art_wnd  ) );
}


void on_sino_advart_trick1_toggled (GtkToggleButton *togglebutton, gpointer user_data)
{
	art_use_smtrick = togglebutton->active;
	dp (2, "Smoothing trick now %d\n",art_use_smtrick);
}

void on_sino_advart_trick2_toggled (GtkToggleButton *togglebutton, gpointer user_data)
{
	art_use_ctrick = togglebutton->active;
	dp (2, "Constraint trick now %d\n",art_use_ctrick);
}

void on_sino_advart_hamming_toggled (GtkToggleButton *togglebutton, gpointer user_data)
{
	art_use_hamming = togglebutton->active;
	dp (2, "Hamming window now %d\n",art_use_hamming);
}

void set_art_refcorr_menu (GtkWidget *widget, gpointer user_data)
{
	art_use_refcorr = *(int*) user_data;
	dp (2, "SART refraction tolerance is now %d\n",art_use_refcorr);
	
}

void art_set_maxiter (GtkWidget *widget, gpointer user_data)
{
GtkAdjustment *adj;

	adj = user_data;
	art_maxiter = adj->value;
	dp (2, "SART maximum iterations %d\n",art_maxiter);
	
}


void art_advanced ()
{
GtkWidget *vbox, *table, *hs1;
GtkWidget *label1, *label2, *label3, *label4, *label5;
GtkWidget *trick1toggle, *trick2toggle, *hamtoggle;
GtkWidget *dismissbutton3;
GtkWidget *art_maxiter; GtkAdjustment *art_maxiter_adj;
GtkWidget *glade_menuitem;
GtkWidget *art_refraction_menu, *art_refraction_menu_menu;
GtkTooltips *tt_smooth, *tt_constraint, *tt_hamming, *tt_refcorr, *tt_maxiter;
int y;


	advanced_art_wnd = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	gtk_object_set_data (GTK_OBJECT (advanced_art_wnd), "art_awnd", advanced_art_wnd);
	gtk_window_set_title (GTK_WINDOW (advanced_art_wnd), "SART Advanced Controls");

	vbox = vbox_in_window (advanced_art_wnd, "artvbox2");
	gtk_container_border_width(GTK_CONTAINER(vbox), 5);
	table = MakeTable (advanced_art_wnd,vbox,"arttab2",2,9);
	y=0;

	label1 = MakeLabel (advanced_art_wnd,table,"Use Smoothing Trick","artlabel1",0,y);
	trick1toggle = Checkbutton_in_table (advanced_art_wnd, table, "      ", "toggle1_tog", 1,y, on_sino_advart_trick1_toggled);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (trick1toggle), art_use_smtrick);
	tt_smooth = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_smooth, trick1toggle, "Check to enable between-iteration smooting.", NULL);
	y++;

	label2 = MakeLabel (advanced_art_wnd,table,"Use Constraint Trick","artlabel2",0,y);
	trick2toggle = Checkbutton_in_table (advanced_art_wnd, table, "      ", "toggle2_tog", 1,y, on_sino_advart_trick2_toggled);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (trick2toggle), art_use_ctrick);
	tt_constraint = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_constraint, trick2toggle, "Check to clamp negative absorbance to 0\n(image value constraints).", NULL);
	y++;

	label3 = MakeLabel (advanced_art_wnd,table,"Use Hamming Window","artlabel3",0,y);
	hamtoggle = Checkbutton_in_table (advanced_art_wnd, table, "      ", "toggle3_tog", 1,y, on_sino_advart_hamming_toggled);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (hamtoggle), art_use_hamming);
	tt_hamming = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_hamming, hamtoggle, "Check to apply Hamming window to SART backprojection.", NULL);
	y++;

	label4 = MakeLabel (advanced_art_wnd,table,"Refraction Tolerance","artlabel4",0,y);
	art_refraction_menu = Make_Optionmenu (advanced_art_wnd, table,"art_refraction_menu",1,y);
	gtk_container_set_border_width (GTK_CONTAINER (art_refraction_menu), 1);
	art_refraction_menu_menu = gtk_menu_new ();
	glade_menuitem = Option_menuitem (art_refraction_menu_menu, GTK_SIGNAL_FUNC(set_art_refcorr_menu), "No correction", (gpointer) &menuopt0);
	glade_menuitem = Option_menuitem (art_refraction_menu_menu, GTK_SIGNAL_FUNC(set_art_refcorr_menu), "Omit refraction zone", (gpointer) &menuopt1);
	glade_menuitem = Option_menuitem (art_refraction_menu_menu, GTK_SIGNAL_FUNC(set_art_refcorr_menu), "Interpolate over refraction zone", (gpointer) &menuopt2);
	gtk_option_menu_set_menu (GTK_OPTION_MENU (art_refraction_menu), art_refraction_menu_menu);
	gtk_option_menu_set_history (GTK_OPTION_MENU (art_refraction_menu), art_use_refcorr);
	tt_refcorr = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_refcorr, art_refraction_menu, "For optical tomography ONLY,\nchoose correction of refraction zones.", NULL);
	y++;

	label5 = MakeLabel (advanced_art_wnd,table,"Maximum iterations","artlabel5",0,y);
	art_maxiter_adj = GTK_ADJUSTMENT (gtk_adjustment_new (10, 1, 25, 1, 10, 0));
	art_maxiter = MakeHScale (art_maxiter_adj,advanced_art_wnd,table,"art_maxiter",1,y);
	gtk_scale_set_digits (GTK_SCALE (art_maxiter), 0);
	tt_maxiter = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_maxiter, art_maxiter, "Limits number of ART/SART iterations.", NULL);


	hs1 = MakeHsep (advanced_art_wnd,vbox, "artsep3",-1,-1);
 	dismissbutton3 = Button_in_table (advanced_art_wnd,vbox,"Dismiss","swbutton3",-1,-1,destroy_advart_wnd);

	gtk_signal_connect (GTK_OBJECT (art_maxiter_adj), "value_changed",
                      GTK_SIGNAL_FUNC (art_set_maxiter), art_maxiter_adj);
	gtk_signal_connect (GTK_OBJECT (advanced_art_wnd), "destroy",
                    	GTK_SIGNAL_FUNC (destroy_advart_wnd),advanced_art_wnd);

	sino_advart_exists = 1;
	gtk_widget_show (advanced_art_wnd);


}



void on_recon_cancelbutton_clicked  (GtkButton *button, gpointer user_data)
{
	if (sino_advart_exists) destroy_advart_wnd (NULL,NULL);

	gtk_widget_destroy (GTK_WIDGET (user_data));
}


void on_recon_okbutton_clicked  (GtkButton *button, gpointer user_data)
{
int res;


	undo_nexus (&mainimg);

	if (sino_realworld)
	{
		if (recon_manualshift)
			ctprepare (&mainimg, sino_logdata, (int)recon_shift_adj->value, 1, recon_pixsize_adj->value, sino_denoise);
		else
			ctprepare (&mainimg, sino_logdata, 0, 0, recon_pixsize_adj->value, sino_denoise);
	}

	res = recon_resolution_adj->value;

	switch (recon_kernel)
	{
		case 1:				/* Backproject, don't filter */
			reconstruct_slice (&mainimg, res, sino_halfview);
			break;
		case 2:				/* Shepp-Logan filter, don't backproject */
			reconstruct_prefilter (&mainimg, res, 1);
			break;
		case 3:				/* Shepp-Logan filtered backprojection */
			reconstruct_prefilter (&mainimg, res, 1);
			reconstruct_slice (&mainimg, res, sino_halfview);
			break;
		case 4:				/* Ramachandran filter, don't backproject */
			reconstruct_prefilter (&mainimg, res, 2);
			break;
		case 5:				/* Ramachandran filtered backprojection */
			reconstruct_prefilter (&mainimg, res, 2);
			reconstruct_slice (&mainimg, res, sino_halfview);
			break;
		case 6:				/* Radon reconstruction using fft */
			fft_backproject (&mainimg, res, sino_halfview, 0);
			break;
		case 7:				/* SART reconstruction */
			if (art_use_hamming)
				art_recon (&mainimg,2, sino_halfview);
			else
				art_recon (&mainimg,0, sino_halfview);
			break;
		case 8:				/* Kaczmarz ART reconstruction */
			art_recon (&mainimg,1, sino_halfview);
			break;
		case 9:				/* Test ART geometry */
			if (sino_halfview)
				wijtest (&mainimg, res, 180/5, 0);
			else
				wijtest (&mainimg, res, 360/5, 0);
			break;
	}
	make_imlib_buffer_from_main_image (cur_slice);
	gtk_widget_destroy (GTK_WIDGET (user_data));
}


void set_recon_kernel_menu (GtkWidget *widget, gpointer user_data)
{
	recon_kernel = *(int*) user_data;
	
}


void set_recon_denoise_menu  (GtkButton *button, gpointer user_data)
{
	sino_denoise =  *(int*) user_data;
	switch (sino_denoise)
	{
		case 0:
			dp (1, "No sinogram denoising\n");
			break;
		case 1:
			dp (1, "Denoising with AMMSE filter\n");
			break;
		case 2:
			dp (1, "Denoising with bilateral filter\n");
			break;
		case 4:
			dp (1, "Denoising with anisotropic diffusion\n");
			break;
		case 3:
			dp (1, "Denoising with wavelet-based BayesShrink\n");
			break;
		default:
			dp (1, "Unknown denoising method\n");
	}

}


void recon_preview (GtkWidget* widget, gpointer user_data)
{
GtkWidget* label4;

	label4 = user_data;
	gtk_label_set_text (GTK_LABEL(label4), "Manual Centering");
	recon_manualshift=1;
	fast_fft_preview (&mainimg, (int)recon_shift_adj->value, sino_realworld, sino_logdata, sino_halfview);

}


void on_sino_artopts_toggled (GtkToggleButton *togglebutton, gpointer user_data)
{
	if (sino_advart_exists>0)
	{
		sino_advart_exists = 0;
		gtk_widget_destroy ( GTK_WIDGET ( advanced_art_wnd ) );
	}
	else
		if (sino_advart_exists==0) art_advanced ();

}



void create_reconstruction_wnd (image_type* img)
{
GtkWidget* wnd;
GtkWidget* vbox;
GtkWidget* table, *hs1, *hs2;
GtkWidget *label1, *label2, *label3, *label4, *label5;
GtkWidget *recon_resolution, *recon_shift, *recon_pixsize;
GtkWidget *fanbeam_checkbutton, *realworld_checkbutton, *logdata_checkbutton, *halfview_checkbutton;
GtkWidget *straighten_checkbutton, *denoise_checkbutton;
GtkWidget *recon_okbutton, *recon_cancelbutton;
GtkWidget *glade_menuitem;
GtkWidget *recon_kernel_menu, *recon_kernel_menu_menu;
GtkWidget *recon_denoise_menu, *recon_denoise_menu_menu;
GtkTooltips *tt_numdets, *tt_recon, *tt_advsart, *tt_fanbeam, *tt_realworld, *tt_islog, *tt_is360,
			*tt_straighten, *tt_center, *tt_pixsize, *tt_sinodenoise;
int y;


	/* Defaults for ART-SART */
	art_use_smtrick=1; art_use_ctrick=1; art_use_hamming=0; art_use_refcorr=0;
	art_maxiter = 10; art_epsmax = 100;
	sino_halfview=0; sino_denoise=0;

	wnd = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	gtk_object_set_data (GTK_OBJECT (wnd), "wnd002", wnd);
	gtk_container_set_border_width (GTK_CONTAINER (wnd), 3);
	gtk_window_set_title (GTK_WINDOW (wnd), "CT Reconstruction");
	gtk_window_set_modal (GTK_WINDOW (wnd), FALSE);
	vbox = vbox_in_window (wnd,"vbox002");
	table = MakeTable (wnd,vbox,"table002",2,12);

	y=0;
	label1 = MakeLabel (wnd,table,"Pixels/Scan","label1",0,y);
	recon_resolution_adj = GTK_ADJUSTMENT (gtk_adjustment_new (img->xmax, 32, 1024, 1, 10, 0));
	recon_resolution = MakeHScale (recon_resolution_adj,wnd,table,"recon_res",1,y);
	gtk_scale_set_digits (GTK_SCALE (recon_resolution), 0);
	tt_numdets = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_numdets, recon_resolution, "Number of detectors (width of sinogram)\nThis slider is obsolete,\nand changing its value may cause unusual effects.", NULL);

	y++;
	label5 = MakeLabel (wnd,table,"Resolution (mm/pixel)","label5",0,y);
	recon_pixsize_adj = GTK_ADJUSTMENT (gtk_adjustment_new (0.5, 0.05, 5.0, 0.05, 0.5, 0));
	recon_pixsize = MakeHScale (recon_pixsize_adj,wnd,table,"recon_pxsiz",1,y);
	gtk_scale_set_digits (GTK_SCALE (recon_pixsize), 2);
	tt_pixsize = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_pixsize, recon_pixsize, "Detector or pixel size in mm.\nThis value scales the sinogram\nto represent effective absorption.\nThis is only used when logarithmizing real-world data.", NULL);

	y++;
	label2 = MakeLabel (wnd,table,"Reconstruction method","label2",0,y);
	recon_kernel_menu = Make_Optionmenu (wnd, table,"recon_kernel_menu",1,y);
	gtk_container_set_border_width (GTK_CONTAINER (recon_kernel_menu), 1);
	recon_kernel_menu_menu = gtk_menu_new ();
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "No reconstruction", (gpointer) &menuopt0);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "Backproject only", (gpointer) &menuopt1);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "Shepp-Logan", (gpointer) &menuopt2);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "Shepp-Logan backproject", (gpointer) &menuopt3);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "Ramachandran", (gpointer) &menuopt4);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "Ramachandran backproject", (gpointer) &menuopt5);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "Radon reconstruction", (gpointer) &menuopt6);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "SART reconstruction", (gpointer) &menuopt7);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "ART reconstruction", (gpointer) &menuopt8);
	glade_menuitem = Option_menuitem (recon_kernel_menu_menu, GTK_SIGNAL_FUNC(set_recon_kernel_menu), "ART geom test", (gpointer) &menuopt9);
	gtk_option_menu_set_menu (GTK_OPTION_MENU (recon_kernel_menu), recon_kernel_menu_menu);
	tt_recon = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_recon, recon_kernel_menu, "Reconstruction method or kernel.\nNote that some methods are for demonstration purposes only.\nChoose a filtered backprojection, Fourier (Radon), or ART/SART reconstruct.\nART and SART can only handle single slices.", NULL);

	y++;
	label2 = MakeLabel (wnd,table,"Advanced SART options","label2a",0,y);
	art_opt_button = Checkbutton_in_table (wnd,table," ","artopt_chkbutton",1,y,on_sino_artopts_toggled);
	tt_advsart = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_advsart, art_opt_button, "Click this button to access the\nadvanced ART/SART options.", NULL);

	y++;
	label3 = MakeLabel (wnd,table,"Fan-beam geometry","label3",0,y);
	fanbeam_checkbutton = Checkbutton_in_table (wnd,table,"(not yet implemented)", "sino_chkbutton",1,y,on_sino_fanbeam_clicked);
	tt_fanbeam = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_fanbeam, fanbeam_checkbutton, "This button has no effect,\nbecause the fan-beam geometry is not implemented.", NULL);

	y++;
	hs1 = MakeHsep (wnd,table, "hsep3aa",0,y); hs1 = MakeHsep (wnd,table, "hsep3ab",1,y);

	y++;
	realworld_checkbutton = Checkbutton_in_table (wnd,table,"Real-world data", "sinorw_chkbutton",0,y,on_sino_realwld_clicked);
	logdata_checkbutton = Checkbutton_in_table (wnd,table,"already logarithmized","sinorw_chkbutton",1,y,on_sino_logdata_clicked);
	tt_realworld = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_realworld, realworld_checkbutton, "When this is checked, data are assumed to be intensity profiles,\nand a conversion to absorbance is performed.", NULL);
	tt_islog = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_islog, logdata_checkbutton, "When this is checked, data are assumed to be log-transformed.\nAbsorbance is computed as I0-I rather than ln(I0/I).", NULL);

	y++;
	halfview_checkbutton = Checkbutton_in_table (wnd,table,"Sinogram is full 360°", "sinofv_chkbutton",0,y,NULL);
	straighten_checkbutton = Checkbutton_in_table (wnd,table,"Straighten sinogram", "sinost_chkbutton",1,y,on_sino_straighten_clicked);
	tt_is360 = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_is360, halfview_checkbutton, "Check this if the sinogram covers a full revolution,\nuncheck if the sinogram was collected over 180 degrees.", NULL);
	tt_straighten = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_straighten, straighten_checkbutton, "Unimplemented function.\nThis button has no effect.", NULL);

	y++;
	recon_denoise_menu = Make_Optionmenu (wnd, table,"recon_denoise_menu",0,y);
	gtk_container_set_border_width (GTK_CONTAINER (recon_denoise_menu), 1);
	recon_denoise_menu_menu = gtk_menu_new ();
	glade_menuitem = Option_menuitem (recon_denoise_menu_menu, GTK_SIGNAL_FUNC(set_recon_denoise_menu), "No denoising", (gpointer) &menuopt0);
	glade_menuitem = Option_menuitem (recon_denoise_menu_menu, GTK_SIGNAL_FUNC(set_recon_denoise_menu), "AMMSE filter", (gpointer) &menuopt1);
	glade_menuitem = Option_menuitem (recon_denoise_menu_menu, GTK_SIGNAL_FUNC(set_recon_denoise_menu), "Bilateral filter", (gpointer) &menuopt2);
	glade_menuitem = Option_menuitem (recon_denoise_menu_menu, GTK_SIGNAL_FUNC(set_recon_denoise_menu), "Wavelet-BayesShrink", (gpointer) &menuopt3);
	gtk_option_menu_set_menu (GTK_OPTION_MENU (recon_denoise_menu), recon_denoise_menu_menu);
	tt_sinodenoise = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_sinodenoise, recon_denoise_menu, "Select a sinogram-domain denoising method", NULL);

	y++;
	label4 = MakeLabel (wnd,table,"Automatic centering","label4",0,y);
	recon_shift_adj = GTK_ADJUSTMENT (gtk_adjustment_new (0, -256, 256, 1, 10, 0));
	recon_shift = MakeHScale (recon_shift_adj,wnd,table,"recon_shift",1,y);
	gtk_scale_set_digits (GTK_SCALE (recon_shift), 0);
	tt_center = gtk_tooltips_new ();
	gtk_tooltips_set_tip (tt_center, recon_shift, "Sinogram centering function.\nLeave unchanged for automatic centering\n(only possible for 360 degree sinograms)\nor change to adjust manually in a preview window.", NULL);

	y++;
	hs2 = MakeHsep (wnd,table, "hsep3ba",0,y); hs2 = MakeHsep (wnd,table, "hsep3bb",1,y);
	y++;
	recon_okbutton = Button_in_table (wnd,table,"OK","recon_okbutton",0,y,NULL);
	recon_cancelbutton = Button_in_table (wnd,table,"Cancel","recon_cancelbutton",1,y,NULL);

	gtk_signal_connect (GTK_OBJECT (recon_shift_adj), "value_changed",
                      GTK_SIGNAL_FUNC (recon_preview), label4);
	gtk_signal_connect (GTK_OBJECT (wnd), "destroy",
                      GTK_SIGNAL_FUNC (on_recon_cancelbutton_clicked), wnd);
	gtk_signal_connect (GTK_OBJECT (recon_cancelbutton), "clicked",
                      GTK_SIGNAL_FUNC (on_recon_cancelbutton_clicked), wnd);
	gtk_signal_connect (GTK_OBJECT (recon_okbutton), "clicked",
                      GTK_SIGNAL_FUNC (on_recon_okbutton_clicked), wnd);
	gtk_signal_connect (GTK_OBJECT (halfview_checkbutton), "toggled",
                      GTK_SIGNAL_FUNC (on_sino_halfview_clicked), halfview_checkbutton);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(halfview_checkbutton), TRUE);		/* Note that this emits a toggle signal and calls the above function */

	gtk_widget_show (wnd);
}




void sino_main ()
{

	art_sintbl=NULL; art_costbl=NULL;
	sino_fanbeam = 0; sino_realworld=0; sino_logdata=0; sino_alg=0; sino_denoise=0;
	create_sinogram_wnd ();
}



void ct_reconstruct_main ()
{
	sino_fanbeam = 0; recon_kernel=0; sino_realworld=0;
	recon_manualshift=0; sino_straighten=0;

	create_reconstruction_wnd (&mainimg);

}



/* Macro interpreter hookup. The syntax is

	imgvar = reconstruct (imgexpr, "method", prep, deltax);
	imgvar = reconstruct (imgexpr, "method", prep, deltax, denoise);
	imgvar is a N x N image, where N is the x-dimension of
	imgexpr. Imgexpr describes a full 360° sinogram. In some cases,
	Method is a string expr which determines the reconstruction method.
	It must be one of "fbp", "radon", "art" or "sart".
		Note: FBP is always with Shepp-Logan kernel
	prep determines sinogram preparation:
		0: No preparation, sinogram is an absorbance sinogram
		1: Sinogram is from a real-world instrument, determine I0 from shoulders and take ln (I/I0)
		2: Sinogram is from a real-world instrument, but is already logarithmized.
				determine I0(log) from shoulders and take I0-I
	deltax is the width of a pixel in mm (floating-point value)
	denoise is an integer that chooses the sinogram-level denoising function:
		0: No denoising
		1: AMMSE
		2: Bilateral filter

*/


variable_type mac_reconstruct (unsigned char* s, int* i)
{
variable_type v,r;
int j, method, halfview, prep, req_denoise;
float deltax;
char fn[20];

	j = *i; 
	v = expr (s,&j, MT_IMG);

	if (s[j++]!=',') macro_error (1,j);
	r = expr (s, &j, MT_STRING);
	strncpy (fn, r.strg, 19);
	free (r.strg);

	/* Default advanced opts - TODO: Need to become accessible */
	
	art_use_refcorr=0;
	art_use_hamming=0;
	art_use_smtrick=1;
	art_use_ctrick=1;
	art_maxiter = 10; art_epsmax = 100;
	halfview = 0;		/* fixed 360 degree sinogram */
	req_denoise = 0;	/* No denoising by default */
	method=0;			/* Silence a compiler warning */


	if (!strncasecmp (fn,"fbp360",6)) 
	{ 
		method=1; halfview=0; 
	}
	else if (!strncasecmp (fn,"fbp180",6)) 
	{ 
		method=1; halfview=1; 
	}
	else if (!strncasecmp (fn,"radon360",8))
	{
		method=2; halfview=0;
	}
	else if (!strncasecmp (fn,"radon180",8))
	{
		method=2; halfview=1;
	}
	else if (!strncasecmp (fn,"art360",6)) 
	{
		method=3; halfview=0;
	}
	else if (!strncasecmp (fn,"art180",6)) 
	{
		method=3; halfview=1;
	}
	else if (!strncasecmp (fn,"sart360",7))
	{
		method=4; halfview=0;
	}
	else if (!strncasecmp (fn,"sart180",7))
	{
		method=4; halfview=1;
	}
	else macro_error (15,j);

	if (s[j++]!=',') macro_error (1,j);
	r = expr (s, &j, MT_INT);
	prep = r.intvalue;

	if (s[j++]!=',') macro_error (1,j);
	r = expr (s, &j, MT_REAL);
	deltax = r.numvalue;

	if (s[j]==',') 			/* Optional denoising parameter  */
	{
		j++;
		r = expr (s,&j,MT_INT);
		req_denoise = r.intvalue;
	}

	while (s[j]!=')') j++;
	*i = j+1;

	if (prep==1)
		ctprepare (&v.imgptr, 0, 0, 1, deltax, req_denoise);
	else if (prep==2)
		ctprepare (&v.imgptr, 1, 0, 1, deltax, req_denoise);



	switch (method)
	{
		case 1:
			reconstruct_prefilter (&v.imgptr, v.imgptr.xmax, 1);
			reconstruct_slice(&v.imgptr, v.imgptr.xmax, halfview);
			break;
		case 2:
			fft_backproject (&v.imgptr, v.imgptr.xmax, halfview, 0);
			break;
		case 3:
			art_recon (&v.imgptr, 0, halfview);
			break;
		case 4:
			art_recon (&v.imgptr, 1, halfview);
			break;
	}

	return v;

}


/* Macro interpreter hookup. The syntax is

	imgvar = sinogram (imgexpr, numexpr);
	imgvar = sinogram (imgexpr, numexpr, strexpr);
	imgexpr is a N x N image, and sinogram creates projections
	with an angular increment described by numexpr.
	Strexpr determines the projection method
	and may be one of "art" or "ray[trace]". By default, 
	this macro function uses the raytrace projection method
	
*/


variable_type mac_sinogram (unsigned char* s, int* i)
{
variable_type v,r;
float anginc;
int j, method;

	j = *i; 
	v = expr (s,&j, MT_IMG);

	if (s[j++]!=',') macro_error (1,j);
	r = expr (s, &j, MT_REAL);		/* Angular increment */
	anginc = r.numvalue;

	method = 1;
	if (s[j]==',') 			/* Optional projection method  */
	{
		j++;
		r = expr (s,&j,MT_STRING);
		if (!strncasecmp (r.strg,"art",3)) method=0;
		else if (!strncasecmp (r.strg,"ray",3)) method=1;
	}

	while (s[j]!=')') j++;
	*i = j+1;

	if (method==0)
		art_project (&v.imgptr, 0, MAX(v.imgptr.xmax,v.imgptr.ymax), anginc);
	else if (method==1)
		make_sinogram_simple (&v.imgptr, 0, MAX(v.imgptr.xmax,v.imgptr.ymax), anginc, 0);
	return v;

}



/**************************************************************

	This is the infrastructure we use to register with the 
	main program. The pointer "menufunction" points to the 
	main entry (the function to be called when the menu is 
	activated). Initmodule is the one function that the 
	main program calls (the name is fixed!) to initialize
	the module and to register the module with the menu.

***************************************************************/

void (*menufunction1)();
void (*menufunction2)();
variable_type (*tomo_recon)(unsigned char* s, int* i);
variable_type (*tomo_sino)(unsigned char* s, int* i);


void module_make_sino()
{

	dp (2, "Tomography module activated\n");
	if (file_ok) sino_main ();
}


void module_reconstruct()
{

	dp (2, "Tomography module activated\n");
	if (file_ok) ct_reconstruct_main ();
}



void initmodule ()
{
GtkWidget* menu;

	dp (1, "Tomography module: Initializer called\n");
	
	/* Register this module with the main menu */
	
	menufunction1 = &module_make_sino;
	menufunction2 = &module_reconstruct;
	menu = dso_menu_register_with_submenu ("Tomography");
	dso_submenu_register (menu, "Make Sinogram", menufunction1);
	dso_submenu_register (menu, "Reconstruction", menufunction2);

	/* Register the module with the macro interpreter */
	
	tomo_recon = &mac_reconstruct;
	register_dso_macro ("RECONSTRUCT",2,MT_IMG,tomo_recon);
	tomo_sino = &mac_sinogram;
	register_dso_macro ("SINOGRAM",2,MT_IMG,tomo_sino);
	

	
}

