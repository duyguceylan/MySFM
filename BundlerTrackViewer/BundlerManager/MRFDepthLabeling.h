/*
 *  MRFDepthLabeling.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 4/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "Image.h"
#include "PerspectiveCamera.h"
#include <string>

using namespace std;

typedef struct NeighborView
{
	int index;
	float viewAngle;
	
	NeighborView(int index_, float viewAngle_) {index = index_; viewAngle = viewAngle_;}
	
	bool operator<(const NeighborView& rhs) const
	{
		return viewAngle > rhs.viewAngle;
	}
};

class MRFDepthLabeling
{
public:
	MRFDepthLabeling(vector<string> &imageFilenames_, vector<PerspectiveCamera> &cameras_);
	void computeDepthMap(int minX, int maxX, int minY, int maxY, vector<int> &pixelIndices, float minDepth_, float maxDepth_, int refIndex, vector<float> &depths);
	
private:
	vector<string> imageFilenames;
	vector<PerspectiveCamera> cameras;
	float minDepth, maxDepth;
	int depthRes;
	static int smoothScale;
	static int currentWidth;
	static int currentHeight;
	static vector<int> labelIndices;
	static vector<vector<int> > neighborRelations;
	static std::vector<float> hCues;
	static std::vector<float> vCues;
	
	static float computeSmoothnessForMRF(int pix1, int pix2, int label1, int label2);
	void findNeighboringViews(int index, vector<NeighborView>& neighboringViews);
	int resizePixelIndices(vector<int> &oldPixelIndices, vector<int> &newPixelIndices, int widthPrev, int heightPrev, int width, int height);
	void resizeLabels(vector<int> &prevVector, vector<bool> &prevVectorFlag, int widthPrev, int heightPrev, int width, int height, int resPrev, int res);
	void resizeEnergy(vector<float> &prevVector, int widthPrev, int heightPrev, int width, int height, int resPrev, int res);
	void computeGradientSensitiveMrfCues(Img &refImg, vector<int> &pixelIndices, int minX, int maxX, int minY, int maxY, vector<float> &vCues, vector<float> &hCues);
	void computeDepthMapForImagePair(int indexRef, int indexNeighbor, int level, Img& refImg, Img& neighborImg, int minX, int maxX, int minY, int maxY,
									 vector<int> &pixelIndices, vector<float> &nccEnergies, int levelRes, bool limitSearch,
									 float depthOffset, vector<int> &prevDepths, vector<bool> &prevDepthAssigned);
};