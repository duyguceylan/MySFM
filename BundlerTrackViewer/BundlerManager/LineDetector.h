/*
 *  LineDetector.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 7/13/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _LINE_DETECTOR_H
#define _LINE_DETECTOR_H

#include "Photo.h"
#include "Line.h"

class LineDetector
{
public:
	LineDetector() {}
	
	void find2DLineSegments(Photo *p, float sigma, float tlow, float thigh, float lengthThreshold, vector<Line> &lineSegments, const char *filename);
	void read2DLineSegments(const char *filename, vector<Line> &lineSegments, float thresh);
	void save2DLineSegments(const char *filename, vector<Line> &lineSegments, float thresh);
	
	//finding edges compatible with user input
	void findLinesIntersectingRegion(vector<Line> &edges2D, vector<Line> &lines, Vec2f center, float regionWidth, float regionHeight);
	int findAvgCompatibleEdge(vector<Line> &lines, Vec2f startPt, Vec2f endPt, Vec3f &avgLine, vector<int> &closeEdges, int boundaryThresh, Vec2f translation = Vec2f(0.0, 0.0));
	int findBestCompatibleEdge(vector<Line> &edges2D, Vec2f startPt, Vec2f endPt, float &score, int boundaryThresh, Vec2f translation = Vec2f(0.0, 0.0));
	
	void findRepeatingPathsWithLineCompatibility(vector<Line> pathLines, vector<Line> &lines, vector<Vec2f> &matches);
private:
	//edge detection
	void cannyEdgeDetection(Photo *p, float sigma, float tlow, float thigh, unsigned char **edge);
	void computeGradientMagnitudes(short int *derX, short int *derY, int width, int height, short int **magnitude);
	void computeDerivatives(float *smoothedImgR, float *smoothedImgG, float *smoothedImgB, int width, int height, short int **derX, short int **derY);
	void followEdges(unsigned char *edgemapptr, short *edgemagptr, short lowval, int width);
	void applyHysteresis(short int *mag, unsigned char *nms, int width, int height, float tlow, float thigh, unsigned char *edge);
	void computeNonMaxSuppresion(short int *magnitude, short int *derX, short int *derY, int width, int height, unsigned char *result);
	
	//edge linking
	void cleanEdgeImage(unsigned char *image, int width, int height);
	void removeIsolatedPoints(unsigned char *edges, int width, int height);
	void linkOpenEdges(unsigned char *edges, int width, int height, int minLength, vector<vector<int> > &xPixels, vector<vector<int> > &yPixels, Img *im);
	void linkClosedEdges(unsigned char *edges, int width, int height, int minLength, vector<vector<int> > &xPixels, vector<vector<int> > &yPixels, Img *im);
	void linkEdgePixels(int width, int height, unsigned char *edges, int minLength, vector<vector<int> > &xPixels, vector<vector<int> > &yPixels);
	
	//line segment computation
	float measureLineDeviation(vector<int> &xPixels, vector<int> &yPixels, int first, int last, int &index);
	void segmentLinkedLists(int width, int height, float lineDevThresh, int minLineLength, vector<vector<int> > &xPixels, vector<vector<int> > &yPixels, vector<vector<Line> > &lineSegments, Img *img);
	
};

#endif