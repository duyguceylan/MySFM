/*
 *  EdgeDetector.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 12/7/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _EDGE_DETECTOR_H
#define _EDGE_DETECTOR_H

#include "Image.h"

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0

class EdgeDetector
{
public:
	EdgeDetector(){}
	
	void cannyEdgeDetection(Img *p, float sigma, float tlow, float thigh, float **smoothedImgR, float **smoothedImgG, float **smoothedImgB, unsigned char **edges);
	void computeGradientMagnitudes(short *derX, short *derY, int width, int height, short **magnitude);
	void computeDerivatives(float *smoothedImgR, float *smoothedImgG, float *smoothedImgB, int width, int height, short **derX, short **derY);
	void followEdges(unsigned char *edgemapptr, short *edgemagptr, short lowval, int width);
	void applyHysteresis(short *mag, unsigned char *nms, int width, int height, float tlow, float thigh, unsigned char *edge);
	void computeNonMaxSuppresion(short *magnitude, short *derX, short *derY, int width, int height, unsigned char *result);
};

#endif