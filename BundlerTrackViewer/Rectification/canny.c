#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define VERBOSE 0
#define BOOSTBLURFACTOR 90.0

typedef double Float; 

void canny(unsigned char *image, int rows, int cols, float sigma,
		   float tlow, float thigh, short int* smoothed, int *edge);
static void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma,
							short int *smoothedim);
static void make_gaussian_kernel(float sigma, Float **kernel, int *windowsize);
static void derrivative_x_y(short int *smoothedim, int rows, int cols,
							short int **delta_x, short int **delta_y);
static void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
						  short int **magnitude);
static void apply_hysteresis(short int *mag, short *gradx, short *grady,unsigned char *nms, int rows, int cols,
							 float tlow, float thigh, int *edge);
static double angle_radians(double x, double y);

static void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
						 unsigned char *result); 
static void suppress_edge_by_curature(const short int*smoothed, int *edgemap, int cols, int rows);
/*******************************************************************************
* PROCEDURE: canny
* PURPOSE: To perform canny edge detection.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void canny(unsigned char *image, int rows, int cols, float sigma,
		   float tlow, float thigh, short int* smoothed, int *edge)
{
	unsigned char *nms;        /* Points that are local maximal magnitude. */
	short int *delta_x,        /* The first derivative image, x-direction. */
		*delta_y,        /* The first derivative image, y-direction. */
		*magnitude;      /* The magnitude of the gadient image.      */


	/****************************************************************************
	* Perform gaussian smoothing on the image using the input standard
	* deviation.
	****************************************************************************/
	if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
	gaussian_smooth(image, rows, cols, sigma, smoothed);

	/****************************************************************************
	* Compute the first derivative in the x and y directions.
	****************************************************************************/
	if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
	derrivative_x_y(smoothed, rows, cols, &delta_x, &delta_y);


	/****************************************************************************
	* Compute the magnitude of the gradient.
	****************************************************************************/
	if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
	magnitude_x_y(delta_x, delta_y, rows, cols, &magnitude);

	/****************************************************************************
	* Perform non-maximal suppression.
	****************************************************************************/
	if(VERBOSE) printf("Doing the non-maximal suppression.\n");
	if((nms = (unsigned char *) calloc(rows*cols,sizeof(unsigned char)))==NULL){
		fprintf(stderr, "Error allocating the nms image.\n");
		exit(1);
	}

	non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);

	/****************************************************************************
	* Use hysteresis to mark the edge pixels.
	****************************************************************************/
	if(VERBOSE) printf("Doing hysteresis thresholding.\n");
	apply_hysteresis(magnitude, delta_x, delta_y, nms, rows, cols, tlow, thigh, edge);
	suppress_edge_by_curature(smoothed, edge, cols, rows);

	/****************************************************************************
	* Free all of the memory that we allocated except for the edge image that
	* is still being used to store out result.
	****************************************************************************/
	free(delta_x);
	free(delta_y);
	free(magnitude);
	free(nms);
}


/*******************************************************************************
* FUNCTION: angle_radians
* PURPOSE: This procedure computes the angle of a vector with components x and
* y. It returns this angle in radians with the answer being in the range
* 0 <= angle <2*PI.
*******************************************************************************/
double angle_radians(double x, double y)
{
	double xu, yu, ang;

#define M_PI  3.14159265358979323846
	xu = fabs(x);
	yu = fabs(y);

	if((xu == 0) && (yu == 0)) return(0);

	ang = atan(yu/xu);

	if(x >= 0){
		if(y >= 0) return(ang);
		else return(2*M_PI - ang);
	}
	else{
		if(y >= 0) return(M_PI - ang);
		else return(M_PI + ang);
	}
}

/*******************************************************************************
* PROCEDURE: magnitude_x_y
* PURPOSE: Compute the magnitude of the gradient. This is the square root of
* the sum of the squared derivative values.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
				   short int **magnitude)
{
	int r, c, pos, sq1, sq2;

	/****************************************************************************
	* Allocate an image to store the magnitude of the gradient.
	****************************************************************************/
	if((*magnitude = (short *) calloc(rows*cols, sizeof(short))) == NULL){
		fprintf(stderr, "Error allocating the magnitude image.\n");
		exit(1);
	}

	for(r=0,pos=0;r<rows;r++){
		for(c=0;c<cols;c++,pos++){
			sq1 = (int)delta_x[pos] * (int)delta_x[pos];
			sq2 = (int)delta_y[pos] * (int)delta_y[pos];
			(*magnitude)[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
		}
	}

}

/*******************************************************************************
* PROCEDURE: derrivative_x_y
* PURPOSE: Compute the first derivative of the image in both the x any y
* directions. The differential filters that are used are:
*
*                                          -1
*         dx =  -1 0 +1     and       dy =  0
*                                          +1
*
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void derrivative_x_y(short int *smoothedim, int rows, int cols,
					 short int **delta_x, short int **delta_y)
{
	int r, c, pos;

	/****************************************************************************
	* Allocate images to store the derivatives.
	****************************************************************************/
	if(((*delta_x) = (short *) calloc(rows*cols, sizeof(short))) == NULL){
		fprintf(stderr, "Error allocating the delta_x image.\n");
		exit(1);
	}
	if(((*delta_y) = (short *) calloc(rows*cols, sizeof(short))) == NULL){
		fprintf(stderr, "Error allocating the delta_x image.\n");
		exit(1);
	}

	/****************************************************************************
	* Compute the x-derivative. Adjust the derivative at the borders to avoid
	* losing pixels.
	****************************************************************************/
	if(VERBOSE) printf("   Computing the X-direction derivative.\n");
	for(r=0;r<rows;r++){
		pos = r * cols;
		(*delta_x)[pos] = smoothedim[pos+1] - smoothedim[pos];
		pos++;
		for(c=1;c<(cols-1);c++,pos++){
			(*delta_x)[pos] = smoothedim[pos+1] - smoothedim[pos-1];
		}
		(*delta_x)[pos] = smoothedim[pos] - smoothedim[pos-1];
	}

	/****************************************************************************
	* Compute the y-derivative. Adjust the derivative at the borders to avoid
	* losing pixels.
	****************************************************************************/
	if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
	for(c=0;c<cols;c++){
		pos = c;
		(*delta_y)[pos] = smoothedim[pos+cols] - smoothedim[pos];
		pos += cols;
		for(r=1;r<(rows-1);r++,pos+=cols){
			(*delta_y)[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
		}
		(*delta_y)[pos] = smoothedim[pos] - smoothedim[pos-cols];
	}
}

/*******************************************************************************
* PROCEDURE: gaussian_smooth
* PURPOSE: Blur an image with a gaussian filter.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma,
					 short int *smoothedim)
{
	int r, c, rr, cc,     /* Counter variables. */
		windowsize,        /* Dimension of the gaussian kernel. */
		center;            /* Half of the windowsize. */
	Float *tempim,        /* Buffer for separable filter gaussian smoothing. */
		*kernel,        /* A one dimensional gaussian kernel. */
		dot,            /* Dot product summing variable. */
		sum;            /* Sum of the kernel weights variable. */

	//FILE* f = fopen("d:\\testn2.txt", "wb");
	/****************************************************************************
	* Create a 1-dimensional gaussian smoothing kernel.
	****************************************************************************/
	if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
	make_gaussian_kernel(sigma, &kernel, &windowsize);
	center = windowsize / 2;


	/****************************************************************************
	* Allocate a temporary buffer image and the smoothed image.
	****************************************************************************/
	if((tempim = (Float *) calloc(rows*cols, sizeof(Float))) == NULL){
		fprintf(stderr, "Error allocating the buffer image.\n");
		exit(1);
	}


	/****************************************************************************
	* Blur in the x - direction.
	****************************************************************************/
	if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
	for(r=0;r<rows;r++){
		for(c=0;c<cols;c++){
			dot = 0;
			sum = 0;
			for(cc=(-center);cc<=center;++cc){
				if(((c+cc) >= 0) && ((c+cc) < cols)){
					sum += kernel[center+cc];
					dot += (((Float)image[r*cols+(c+cc)]) * kernel[center+cc]);
				}
			}
			tempim[r*cols+c] = (dot/sum);

			//printf("%f,%f, %f\n", dot, sum, tempim[r*cols+c]);
		}
	}

	//fclose(f);	exit(0);
	/****************************************************************************
	* Blur in the y - direction.
	****************************************************************************/
	if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
	for(c=0;c<cols;c++){
		for(r=0;r<rows;r++){
			sum = 0.0;
			dot = 0.0;
			for(rr=(-center);rr<=center;rr++){
				if(((r+rr) >= 0) && ((r+rr) < rows)){
					dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
					sum += kernel[center+rr];
				}
			}
			smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
			/*smoothedim[r*cols+c] = (short int)(BOOSTBLURFACTOR * tempim[r * cols + c] + 0.5);*/
		}
	}

	free(tempim);
	free(kernel);
}

/*******************************************************************************
* PROCEDURE: make_gaussian_kernel
* PURPOSE: Create a one dimensional gaussian kernel.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void make_gaussian_kernel(float sigma, Float **kernel, int *windowsize)
{
	int i, center;
	Float x, fx, sum=0.0;

	*windowsize = 1 + 2 * ((int)ceil(2.5 * sigma));
	center = (*windowsize) / 2;

	if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);
	if((*kernel = (Float *) calloc((*windowsize), sizeof(Float))) == NULL){
		fprintf(stderr, "Error callocing the gaussian kernel array.\n");
		exit(1);
	}

	for(i=0;i<(*windowsize);i++){
		x = (Float)(i - center);
		fx = (Float)( pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853)));
		(*kernel)[i] = fx;
		sum += fx;
	}

	for(i=0;i<(*windowsize);i++) (*kernel)[i] /= sum;

	if(VERBOSE){
		printf("The filter coefficients are:\n");
		for(i=0;i<(*windowsize);i++)
			printf("kernel[%d] = %f\n", i, (*kernel)[i]);
	}
}


/*******************************************************************************
* FILE: hysteresis.c
* This code was re-written by Mike Heath from original code obtained indirectly
* from Michigan State University. heath@csee.usf.edu (Re-written in 1996).
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#define VERBOSE 0

#define NOEDGE 0
#define POSSIBLE_EDGE 1
#define EDGE_BREAK 2
#define EDGE 3

/*******************************************************************************
* PROCEDURE: follow_edges
* PURPOSE: This procedure edges is a recursive routine that traces edgs along
* all paths whose magnitude values remain above some specifyable lower
* threshhold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void follow_edges(int *edgemapptr, short int *edgemagptr, int lowval, int cols)
{
	short int *tempmagptr;
	int *tempmapptr;
	int i;
	int x[8] = {1,1,0,-1,-1,-1,0,1},
		y[8] = {0,1,1,1,0,-1,-1,-1};

	for(i=0;i<8;i++)
	{
		tempmapptr = edgemapptr - y[i]*cols + x[i];
		tempmagptr = edgemagptr - y[i]*cols + x[i];

		if( (*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval))
		{
			*tempmapptr =  *edgemapptr;
			follow_edges(tempmapptr,tempmagptr, lowval, cols);
		}
	}
}


int follow_straight_edges(int* edgemap_ptr, short int* mag_ptr,
						   short int* gx_ptr, short int* gy_ptr, int lowval, int cols)
{
	short int *mag, * gx, * gy;
	int *map;
	int i, offset, count = 0;
	char x[8] = {1,1,0,-1,-1,-1,0,1},
		y[8] = {0,1,1,1,0,-1,-1,-1};
	char tt[8] = {0, 0, 0, 0, 0, 0, 0, 0};

	for(i=0;i<8;i++)
	{
		offset = - y[i]*cols + x[i];
		map = edgemap_ptr + offset;
		mag = mag_ptr + offset;
		gx = gx_ptr + offset;
		gy = gy_ptr + offset;
		if( (*map == POSSIBLE_EDGE) && (*mag > lowval))
		{
			if(4 * abs(gx_ptr[0] * gy[0] - gy_ptr[0] * gx[0]) <  mag[0] * mag_ptr[0])
			{
				*map =  *edgemap_ptr;
				count += follow_straight_edges(map,mag, gx, gy, lowval, cols);
			}else
			{
				tt[i] = 1;
			}
		}
	}
	for(i = 0; i < 8; i++)
	{
		if(tt[i])
		{
			offset = - y[i]*cols + x[i];
			map = edgemap_ptr + offset;
			mag = mag_ptr + offset;
			gx = gx_ptr + offset;
			gy = gy_ptr + offset;
			count ++;
			*map = *edgemap_ptr + count;
			count += follow_straight_edges(map,mag, gx, gy, lowval, cols) ;
		}
	}
	return count;
}
///
void suppress_edge_by_curature(const short int*smoothed, int *edgemap, int cols, int rows)
{
	int i, j;
	const short* smoothed_ptr;
	int* edgemap_ptr;
	int dxx, dxy, dyy, va, vb, filtered  = 0;
	for(i = 1; i < rows -1; ++i)
	{
		edgemap_ptr = edgemap + i * cols + 1;
		smoothed_ptr = smoothed + i * cols + 1;
		for(j = 1; j < cols -1; ++j, edgemap_ptr++, smoothed_ptr++)
		{
			if(*edgemap_ptr < EDGE) continue;
			dxx = 4 *( smoothed_ptr[1] + smoothed_ptr[-1] - 2 * smoothed_ptr[0]);
			dyy = 4 *( smoothed_ptr[cols] + smoothed_ptr[-cols] - 2 * smoothed_ptr[0]);
			dxy = (smoothed_ptr[cols + 1] + smoothed_ptr[-cols -1] - smoothed_ptr[cols -1] - smoothed_ptr[1-cols]);
			va = (dxx + dyy) * (dxx + dyy);
			vb =(dxx * dyy - dxy * dxy);
			if(vb > 0 && va < 4 * vb) 
			{
				*edgemap_ptr = NOEDGE;
                filtered ++;
			}
		}
	}
}

/*******************************************************************************
* PROCEDURE: apply_hysteresis
* PURPOSE: This routine finds edges that are above some high threshhold or
* are connected to a high pixel by a path of pixels greater than a low
* threshold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void apply_hysteresis(short int *mag, short *gradx, short *grady, unsigned char *nms, int rows, int cols,
					  float tlow, float thigh, int *edge)
{
	int r, c, pos, numedges, highcount, edge_count;
	int lowthreshold, highthreshold, hist[32768];
	short int maximum_mag;

	/****************************************************************************
	* Initialize the edge map to possible edges everywhere the non-maximal
	* suppression suggested there could be an edge except for the border. At
	* the border we say there can not be an edge because it makes the
	* follow_edges algorithm more efficient to not worry about tracking an
	* edge off the side of the image.
	****************************************************************************/
	for(r=0,pos=0;r<rows;r++){
		for(c=0;c<cols;c++,pos++){
			if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
			else edge[pos] = NOEDGE;
		}
	}

	for(r=0,pos=0;r<rows;r++,pos+=cols){
		edge[pos] = NOEDGE;
		edge[pos+cols-1] = NOEDGE;
	}
	pos = (rows-1) * cols;
	for(c=0;c<cols;c++,pos++){
		edge[c] = NOEDGE;
		edge[pos] = NOEDGE;
	}

	/****************************************************************************
	* Compute the histogram of the magnitude image. Then use the histogram to
	* compute hysteresis thresholds.
	****************************************************************************/
	for(r=0;r<32768;r++) hist[r] = 0;
	for(r=0,pos=0;r<rows;r++){
		for(c=0;c<cols;c++,pos++){
			if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
		}
	}

	/****************************************************************************
	* Compute the number of pixels that passed the nonmaximal suppression.
	****************************************************************************/
	for(r=1,numedges=0;r<32768;r++){
		if(hist[r] != 0) maximum_mag = r;
		numedges += hist[r];
	}

	highcount = (int)(numedges * thigh + 0.5);

	/****************************************************************************
	* Compute the high threshold value as the (100 * thigh) percentage point
	* in the magnitude of the gradient histogram of all the pixels that passes
	* non-maximal suppression. Then calculate the low threshold as a fraction
	* of the computed high threshold value. John Canny said in his paper
	* "A Computational Approach to Edge Detection" that "The ratio of the
	* high to low threshold in the implementation is in the range two or three
	* to one." That means that in terms of this implementation, we should
	* choose tlow ~= 0.5 or 0.33333.
	****************************************************************************/
	r = 1;
	numedges = hist[1];
	while((r<(maximum_mag-1)) && (numedges < highcount)){
		r++;
		numedges += hist[r];
	}
	highthreshold = r;
	lowthreshold = (int)(highthreshold * tlow + 0.5);

	if(VERBOSE){
		printf("The input low and high fractions of %f and %f computed to\n",
			tlow, thigh);
		printf("magnitude of the gradient threshold values of: %d %d\n",
			lowthreshold, highthreshold);
	}

	/****************************************************************************
	* This loop looks for pixels above the highthreshold to locate edges and
	* then calls follow_edges to continue the edge.
	****************************************************************************/
	edge_count = 0;
	for(r=0,pos=0;r<rows;r++){
		for(c=0;c<cols;c++,pos++){
			if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
				edge[pos] = EDGE + edge_count;
				//follow_edges((edge+pos), (mag+pos), lowthreshold, cols);edge_count ++;
				edge_count += (1 + follow_straight_edges(edge + pos, mag + pos, gradx + pos, grady + pos, lowthreshold, cols));

			}
		}
	}

	/****************************************************************************
	* Set all the remaining possible edges to non-edges.
	****************************************************************************/
	for(r=0,pos=0;r<rows;r++){
		for(c=0;c<cols;c++,pos++) if(edge[pos] < EDGE) edge[pos] = NOEDGE;
	}
}

/*******************************************************************************
* PROCEDURE: non_max_supp
* PURPOSE: This routine applies non-maximal suppression to the magnitude of
* the gradient image.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
				  unsigned char *result) 
{
	int rowcount, colcount,count;
	short *magrowptr,*magptr;
	short *gxrowptr,*gxptr;
	short *gyrowptr,*gyptr,z1,z2;
	short m00,gx,gy;
	float mag1,mag2,xperp,yperp;
	unsigned char *resultrowptr, *resultptr;


	/****************************************************************************
	* Zero the edges of the result image.
	****************************************************************************/
	for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1); 
		count<ncols; resultptr++,resultrowptr++,count++){
			*resultrowptr = *resultptr = (unsigned char) 0;
	}

	for(count=0,resultptr=result,resultrowptr=result+ncols-1;
		count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
			*resultptr = *resultrowptr = (unsigned char) 0;
	}

	/****************************************************************************
	* Suppress non-maximum points.
	****************************************************************************/
	for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
		gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
		rowcount<nrows-2; 
	rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
		resultrowptr+=ncols)
	{   
		for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
			resultptr=resultrowptr;colcount<ncols-2; 
			colcount++,magptr++,gxptr++,gyptr++,resultptr++)
		{   
			m00 = *magptr;
			if(m00 == 0){
				*resultptr = (unsigned char) NOEDGE;
				continue;
			}
			else{
				xperp = -(gx = *gxptr)/((float)m00);
				yperp = (gy = *gyptr)/((float)m00);
			}

			//compute curfature


			if(gx >= 0)
			{
				if(gy >= 0)
				 {
					 if (gx >= gy)
					 {  
						 /* 111 */
						 /* Left point */
						 z1 = *(magptr - 1);
						 z2 = *(magptr - ncols - 1);

						 mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

						 /* Right point */
						 z1 = *(magptr + 1);
						 z2 = *(magptr + ncols + 1);

						 mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
					 }
					 else
					 {    
						 /* 110 */
						 /* Left point */
						 z1 = *(magptr - ncols);
						 z2 = *(magptr - ncols - 1);

						 mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

						 /* Right point */
						 z1 = *(magptr + ncols);
						 z2 = *(magptr + ncols + 1);

						 mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp; 
					 }
				 }
				 else
				 {
					 if (gx >= -gy)
					 {
						 /* 101 */
						 /* Left point */
						 z1 = *(magptr - 1);
						 z2 = *(magptr + ncols - 1);

						 mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

						 /* Right point */
						 z1 = *(magptr + 1);
						 z2 = *(magptr - ncols + 1);

						 mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
					 }
					 else
					 {    
						 /* 100 */
						 /* Left point */
						 z1 = *(magptr + ncols);
						 z2 = *(magptr + ncols - 1);

						 mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

						 /* Right point */
						 z1 = *(magptr - ncols);
						 z2 = *(magptr - ncols + 1);

						 mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp; 
					 }
				 }
			}
			else
			{
				if ((gy = *gyptr) >= 0)
				{
					if (-gx >= gy)
					{          
						/* 011 */
						/* Left point */
						z1 = *(magptr + 1);
						z2 = *(magptr - ncols + 1);

						mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

						/* Right point */
						z1 = *(magptr - 1);
						z2 = *(magptr + ncols - 1);

						mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
					}
					else
					{
						/* 010 */
						/* Left point */
						z1 = *(magptr - ncols);
						z2 = *(magptr - ncols + 1);

						mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

						/* Right point */
						z1 = *(magptr + ncols);
						z2 = *(magptr + ncols - 1);

						mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
					}
				}
				else
				{
					if (-gx > -gy)
					{
						/* 001 */
						/* Left point */
						z1 = *(magptr + 1);
						z2 = *(magptr + ncols + 1);

						mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

						/* Right point */
						z1 = *(magptr - 1);
						z2 = *(magptr - ncols - 1);

						mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
					}
					else
					{
						/* 000 */
						/* Left point */
						z1 = *(magptr + ncols);
						z2 = *(magptr + ncols + 1);

						mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

						/* Right point */
						z1 = *(magptr - ncols);
						z2 = *(magptr - ncols - 1);

						mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
					}
				}
			} 

			/* Now determine if the current point is a maximum point */

			if ((mag1 > 0.0) || (mag2 > 0.0))
			{
				*resultptr = (unsigned char) NOEDGE;
			}
			else
			{    
				if (mag2 == 0.0)
					*resultptr = (unsigned char) NOEDGE;
				else
					*resultptr = (unsigned char) POSSIBLE_EDGE;
			}
		} 
	}
}
