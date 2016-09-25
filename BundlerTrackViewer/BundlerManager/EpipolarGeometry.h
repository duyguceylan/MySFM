/*
 *  EpipolarGeometry.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _EPIPOLAR_GEOMETRY_H
#define _EPIPOLAR_GEOMETRY_H

#include <vector>
#include "Common.h"
#include "Plane.h"
#include "PerspectiveCamera.h"
#include "Image.h"

class EpipolarGeometry
{
private:
	static void fmatrixResiduals(const int *m, const int *n, const double *x, double *fvec, int *iflag);
	
	void refineFmatrixNonlinearMatches(int numPts, std::vector<Vec3f> &rightPoints, std::vector<Vec3f> &leftPoints, Matrix3f &F, Matrix3f &FOut);
	
	int estimateFmatrixLinear(std::vector<Vec3f> &rightPoints, std::vector<Vec3f> &leftPoints, Matrix3f &Fout);
	
	int estimateFmatrixRansacMatches(int numPts, vector<Vec3f> &aPts, vector<Vec3f> &bPts, int numTrials, double threshold, 
									 double successRatio, Matrix3f &F, vector<int> &inliers);
	
public:
	EpipolarGeometry() {}
	
	//fundemental matrix
	static float fmatrixComputeResidual(Matrix3f F, Vec3f l, Vec3f r);
	
	std::vector<int> estimateFMatrix(std::vector<keypt_t> &k1, std::vector<keypt_t> &k2, std::vector<KeypointMatch> matches, 
									 int num_trials, double threshold, Matrix3f &F);
	
	std::vector<int> estimateFMatrix(int noPts, std::vector<Vec3f> &k1, std::vector<Vec3f> &k2,
									 int num_trials, double threshold, Matrix3f &F);
	
	bool areFMatricesEqual(Matrix3f F1, Matrix3f F2, int width, int height, vector<Vec2f> &pts, float threshold);
	
	//homography matrix
	bool estimateHomography(int noPts, std::vector<Vec2f> &k1, std::vector<Vec2f> &k2, Matrix3f &H);
	
	//rectification
	void estimateRectifiedTransformed(int numPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts, float &scale, Vec2f &translation);
	
	int estimateRectifiedTransformedRansacMatches(int numPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts, int numTrials, double threshold, 
												  double successRatio, float &scale, Vec2f &translation, vector<int> &inliers);
	
	void estimateRectifiedTransformedRansacMatches(int numPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts, int numTrials, double threshold,
												   vector<float> &scale, vector<Vec2f> &translation, vector<vector<int> > &inliers);
	
	//pose estimation
	int estimatePose(int noPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts, Matrix3f K1, Matrix3f K2, Matrix3f &R, Vec3f &t, vector<int> &inliers);
	
	bool findRotationAndTranslationFromEssentialMatrix(Matrix3f E, Matrix3f K, Vec2f p1, Vec2f p2, Matrix3f &R, Vec3f &t);
	Vec3f triangulate(Vec2f p, Vec2f q, Matrix3f K, Matrix3f R0, Vec3f t0, Matrix3f R1, Vec3f t1, float &error); 
	
	bool isIdentityRotation(Matrix3f R);
	
	//triangulation
	Vec3f computeClosest3DPoint(Vec2f refPixel, Vec2f neighborPixel, PerspectiveCamera &refCam, PerspectiveCamera &neighborCam, int level);
	void findSamplesOnEpipolarLine(Vec3f pixelPos, int pixelIndex, PerspectiveCamera refCam, PerspectiveCamera neighborCam, Img& refImg, Img& neighborImg, int level,
								   float minDepth, float maxDepth, int depthRes, vector<float> &nccEnergies,int minLabelIndex, int totalLabelNo);	
private:
	static float globalScale;
	static std::vector<Vec3f> globalIns;
	static std::vector<Vec3f> globalOuts;
};

#endif