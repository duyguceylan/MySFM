/*
 *  EdgeDetector.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 12/7/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "EdgeDetector.h"


void EdgeDetector::cannyEdgeDetection(Img *p, float sigma, float tlow, float thigh, 
									  float **smoothedImgR, float **smoothedImgG, float **smoothedImgB,
									  unsigned char **edges)
{
	int width = p->width();
	int height = p->height();
	
	unsigned char *nms;
	short *derX;
	short *derY;
	short *magnitude;
	
	(*edges) = new unsigned char[width*height];
	nms = new unsigned char[width*height];
	
	//perform gaussian smoothing
	p->gaussianSmooth(sigma, smoothedImgR, smoothedImgG, smoothedImgB);
	
	//compute gradients
	computeDerivatives(*smoothedImgR, *smoothedImgG, *smoothedImgB, width, height, &derX, &derY);
	
	//compute gradient magnitudes
	computeGradientMagnitudes(derX, derY, width, height, &magnitude);
	
	//compute non-max-suppresion
	computeNonMaxSuppresion(magnitude, derX, derY, width, height, nms);
	
	//apply hysteresis
	applyHysteresis(magnitude, nms, width, height, tlow, thigh, *edges);
	
	delete [] derX;
	delete [] derY;
	delete [] magnitude;
	delete [] nms;	
}

void EdgeDetector::computeGradientMagnitudes(short *derX, short *derY, int width, int height, short **magnitude)
{
	int x, y, pos, sq1, sq2;
	
	(*magnitude) = new short [width*height];
	
	for(y=0,pos=0;y<height;y++)
	{
		for(x=0;x<width;x++,pos++)
		{
			sq1 = (int)derX[pos] * (int)derX[pos];
			sq2 = (int)derY[pos] * (int)derY[pos];
			(*magnitude)[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
		}
	}
	
}

void EdgeDetector::computeDerivatives(float *smoothedImgR, float *smoothedImgG, float *smoothedImgB, int width, int height, short **derX, short **derY)
{
	int x, y, pos;
	
	(*derX) = new short [width*height];
	(*derY) = new short [width*height];
	
	float a, b, c;
	
	//x-derivative
	for(y=0;y<height;y++)
	{
		pos = y * width;
		a = smoothedImgR[pos+1] - smoothedImgR[pos];
		b = smoothedImgG[pos+1] - smoothedImgG[pos];
		c = smoothedImgB[pos+1] - smoothedImgB[pos];
		if(a>b && a>c)
			(*derX)[pos] = a*255;
		else if(b>a && c>a)
			(*derX)[pos] = b*255;
		else
			(*derX)[pos] = c*255;
		//(*derX)[pos] = smoothedImg[pos+1] - smoothedImg[pos]) * 255;
		pos++;
		for(x=1;x<(width-1);x++,pos++)
		{
			a = smoothedImgR[pos+1] - smoothedImgR[pos-1];
			b = smoothedImgG[pos+1] - smoothedImgG[pos-1];
			c = smoothedImgB[pos+1] - smoothedImgB[pos-1];
			if(a>b && a>c)
				(*derX)[pos] = a*255;
			else if(b>a && c>a)
				(*derX)[pos] = b*255;
			else
				(*derX)[pos] = c*255;
			//(*derX)[pos] = (smoothedImg[pos+1] - smoothedImg[pos-1]) * 255;
		}
		a = smoothedImgR[pos] - smoothedImgR[pos-1];
		b = smoothedImgG[pos] - smoothedImgG[pos-1];
		c = smoothedImgB[pos] - smoothedImgB[pos-1];
		if(a>b && a>c)
			(*derX)[pos] = a*255;
		else if(b>a && c>a)
			(*derX)[pos] = b*255;
		else
			(*derX)[pos] = c*255;
		//(*derX)[pos] = (smoothedImg[pos] - smoothedImg[pos-1]) * 255;
	}
	
	//y-derivative
	for(x=0;x<width;x++)
	{
		pos = x;
		a = smoothedImgR[pos+width] - smoothedImgR[pos];
		b = smoothedImgG[pos+width] - smoothedImgG[pos];
		c = smoothedImgB[pos+width] - smoothedImgB[pos];
		if(a>b && a>c)
			(*derY)[pos] = a*255;
		else if(b>a && c>a)
			(*derY)[pos] = b*255;
		else
			(*derY)[pos] = c*255;
		//(*derY)[pos] = (smoothedImg[pos+width] - smoothedImg[pos])*255;
		pos += width;
		for(y=1;y<(height-1);y++,pos+=width)
		{
			a = smoothedImgR[pos+width] - smoothedImgR[pos-width];
			b = smoothedImgG[pos+width] - smoothedImgG[pos-width];
			c = smoothedImgB[pos+width] - smoothedImgB[pos-width];
			if(a>b && a>c)
				(*derY)[pos] = a*255;
			else if(b>a && c>a)
				(*derY)[pos] = b*255;
			else
				(*derY)[pos] = c*255;
			//(*derY)[pos] = (smoothedImg[pos+width] - smoothedImg[pos-width])*255;
		}
		a = smoothedImgR[pos] - smoothedImgR[pos-width];
		b = smoothedImgG[pos] - smoothedImgG[pos-width];
		c = smoothedImgB[pos] - smoothedImgB[pos-width];
		if(a>b && a>c)
			(*derY)[pos] = a*255;
		else if(b>a && c>a)
			(*derY)[pos] = b*255;
		else
			(*derY)[pos] = c*255;
		//(*derY)[pos] = (smoothedImg[pos] - smoothedImg[pos-width])*255;
	}
}

void EdgeDetector::followEdges(unsigned char *edgemapptr, short *edgemagptr, short lowval, int width)
{
	short *tempmagptr;
	unsigned char *tempmapptr;
	int i;
	int x[8] = {1,1,0,-1,-1,-1,0,1}, y[8] = {0,1,1,1,0,-1,-1,-1};
	
	for(i=0;i<8;i++){
		tempmapptr = edgemapptr - y[i]*width + x[i];
		tempmagptr = edgemagptr - y[i]*width + x[i];
		
		if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
			*tempmapptr = (unsigned char) EDGE;
			followEdges(tempmapptr,tempmagptr, lowval, width);
		}
	}
}

void EdgeDetector::applyHysteresis(short *mag, unsigned char *nms, int width, int height, float tlow, float thigh, unsigned char *edge)
{
	int x, y, pos, numedges, highcount, lowthreshold, highthreshold;
	int hist[32768];
	short maximum_mag;
	
	for(y=0,pos=0;y<height;y++)
	{
		for(x=0;x<width;x++,pos++)
		{
			if(nms[pos] == POSSIBLE_EDGE) 
				edge[pos] = POSSIBLE_EDGE;
			else edge[pos] = NOEDGE;
		}
	}
	
	//boundary
	for(y=0,pos=0;y<height;y++,pos+=width)
	{
		edge[pos] = NOEDGE;
		edge[pos+width-1] = NOEDGE;
	}
	pos = (height-1) * width;
	for(x=0;x<width;x++,pos++){
		edge[x] = NOEDGE;
		edge[pos] = NOEDGE;
	}
	
	//histogram
	for(int r=0;r<32768;r++) hist[r] = 0;
	for(y=0,pos=0;y<height;y++)
	{
		for(x=0;x<width;x++,pos++)
		{
			if(edge[pos] == POSSIBLE_EDGE) 
				hist[mag[pos]]++;
		}
	}
	
	numedges = 0;
	for(int r=1;r<32768;r++)
	{
		if(hist[r] != 0) 
			maximum_mag = r;
		numedges += hist[r];
	}
	
	highcount = (int)(numedges * thigh + 0.5);
	
	int r = 1;
	numedges = hist[1];
	while((r<(maximum_mag-1)) && (numedges < highcount)){
		r++;
		numedges += hist[r];
	}
	highthreshold = r;
	lowthreshold = (int)(highthreshold * tlow + 0.5);
	
	for(y=0,pos=0;y<height;y++){
		for(x=0;x<width;x++,pos++){
			if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
				edge[pos] = EDGE;
				followEdges((edge+pos), (mag+pos), lowthreshold, width);
			}
		}
	}
}

void EdgeDetector::computeNonMaxSuppresion(short *magnitude, short *derX, short *derY, int width, int height, unsigned char *result)
{
    int rowcount, colcount,count;
    short *magrowptr,*magptr;
    short *gxrowptr,*gxptr;
    short *gyrowptr,*gyptr,z1,z2;
    short m00,gx,gy;
    float mag1,mag2,xperp,yperp;
    unsigned char *resultrowptr, *resultptr;
	
	//Zero the first and last rows of the result image.
	for(count=0,resultrowptr=result,resultptr=result+width*(height-1); 
        count<width; resultptr++,resultrowptr++,count++)
	{
        *resultrowptr = *resultptr = (unsigned char) 0;
    }
	
	//Zero the first and last columns of the result image
    for(count=0,resultptr=result,resultrowptr=result+width-1;
        count<height; count++,resultptr+=width,resultrowptr+=width)
	{
        *resultptr = *resultrowptr = (unsigned char) 0;
    }
	
	//Suppress non-maximum points.
	for(rowcount=1,magrowptr=magnitude+width+1,gxrowptr=derX+width+1,gyrowptr=derY+width+1,resultrowptr=result+width+1;
		rowcount<height-2; 
		rowcount++,magrowptr+=width,gyrowptr+=width,gxrowptr+=width,resultrowptr+=width)
	{   
		for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,resultptr=resultrowptr;
			colcount<width-2; 
			colcount++,magptr++,gxptr++,gyptr++,resultptr++)
		{   
			m00 = *magptr;
			if(m00 == 0)
			{
				*resultptr = (unsigned char) NOEDGE;
			}
			else
			{
				xperp = -(gx = *gxptr)/((float)m00);
				yperp = (gy = *gyptr)/((float)m00);
			}
			
			if(gx >= 0)
			{
				if(gy >= 0)
				{
                    if (gx >= gy)
                    {  
                        /* 111 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - width - 1);
						
                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                        
                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + width + 1);
						
                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {    
                        /* 110 */
                        /* Left point */
                        z1 = *(magptr - width);
                        z2 = *(magptr - width - 1);
						
                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;
						
                        /* Right point */
                        z1 = *(magptr + width);
                        z2 = *(magptr + width + 1);
						
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
                        z2 = *(magptr + width - 1);
						
                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;
						
                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - width + 1);
						
                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {    
                        /* 100 */
                        /* Left point */
                        z1 = *(magptr + width);
                        z2 = *(magptr + width - 1);
						
                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;
						
                        /* Right point */
                        z1 = *(magptr - width);
                        z2 = *(magptr - width + 1);
						
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
                        z2 = *(magptr - width + 1);
						
                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;
						
                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + width - 1);
						
                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /* 010 */
                        /* Left point */
                        z1 = *(magptr - width);
                        z2 = *(magptr - width + 1);
						
                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;
						
                        /* Right point */
                        z1 = *(magptr + width);
                        z2 = *(magptr + width - 1);
						
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
                        z2 = *(magptr + width + 1);
						
                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;
						
                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - width - 1);
						
                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /* 000 */
                        /* Left point */
                        z1 = *(magptr + width);
                        z2 = *(magptr + width + 1);
						
                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;
						
                        /* Right point */
                        z1 = *(magptr - width);
                        z2 = *(magptr - width - 1);
						
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
				{
                    *resultptr = (unsigned char) NOEDGE;
				}
                else
				{
                    *resultptr = (unsigned char) POSSIBLE_EDGE;
				}
            }
        } 
    }
}
