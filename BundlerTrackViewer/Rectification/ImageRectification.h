/*
 *  ImageRectification.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 12/6/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _IMAGE_RECTIFICATION_H
#define _IMAGE_RECTIFICATION_H

#include <vector>
#include "Image.h"
#include "MyMatrix.h"

#define ORTBINS 36

struct Line2D
{
	float a, b, c;
	Line2D(float xa, float xb, float xc) {a = xa; b = xb; c = xc;}
	Line2D(const float* line) {a = line[0]; b = line[1]; c = line[2];}
	Line2D(const double* line) {a = (float)line[0]; b = (float)line[1]; c = (float)line[2];}
	void SetLine2DU(const double pt[3]) 
	{
		double s = sqrt(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]);
		//double s = 1.0; 
		a = float(pt[0] / s); b = float(pt[1] / s); c = float(pt[2] / s);
	}
	operator float* (){return &a;}
	double DistVP(double x, double y) 
	{
		x = (a - c * x);	y = (b - c * y);
		return sqrt(x * x + y * y) / c;
	}
	double GetAngleCOS(double vp[3])
	{
		return (vp[0] * a + vp[1] * b) /sqrt((vp[0] * vp[0] + vp[1] * vp[1]) * 	(a * a + b * b));
	}
	void GetLine2D(double pt[])
	{
		pt[0] = a; pt[1] = b; pt[2] = c;
	}
	double CompareLine2D(Line2D& l)
	{
		return fabs(a * l.a + b * l.b + c * l.c);
	}
	Line2D(){}
	~Line2D(){}
};

class ImageRectification
{
public:
	ImageRectification ();
	
	void setImage(Img *_im) {im = _im; }
	Matrix3f rectifyImage(Img **data, Vec2f &vanishingPtDist, bool &secondHorizontalVPFound);
	Matrix3f rectifyImageWrtSecondHorizontalDirection(Img **data, Vec2f &vanishingPtDist);
	
	void computeEdgeMap();
	void processEdgeMap();
	void computeDominantOrientations();
	void computeVanishingPoints(Vec2f &vanishingPtDist);
	int computeVanishingPoint(int num, double line[][3], double vp[], int mask[], double th);
	void findSecondHorizontalVP();
	void setHomographyFromVP();
	void computeRectifiedImage(float magnification, Img **data);
	void generateRectifiedImage(int newWidth, int newHeight, Img *data);
	void getRectifiedBoundingBox(int offsetX[2], int offsetY[2], Line2D* vps);
	
	void adjustRectificationMatrix(double dx, double dy);
	void transformPointToOriginal(float x, float y, float &tx, float &ty);
	int computeLineIntersection(double line1[3], double line2[3], double pt[3]);
	int linesDuplicate(double a[3], double b[3]);
	
private:
	Img *im;
	int width, height;
	
	Matrix3f H, Hinv;
	
	unsigned char* imageData;
	short int* smoothData;
	/*unsigned char*/int* edgeData;
	float *smoothedR;
	float *smoothedG;
	float *smoothedB;
	
	vector<Line2D> vPoints;
	
	vector<int> edgePoints;
	vector<Line2D> edgeLines;
	vector<float> edgeOrientations;
	vector<float> edgeMagnitudes;
	
	vector<float> orientationHist;
	vector<int> orientationDom;
	vector<int> orientationLineCount;
	
	bool secondHorizontalDirection;
};

#endif