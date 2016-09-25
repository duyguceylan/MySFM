/*
 *  TextureManager.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 4/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TEXTURE_MANAGER_H
#define TEXTURE_MANAGER_H

#include <vector>
#include "Color.h"

using namespace std;

class TextureManager
{
public:
	static vector<vector<Color> > colors;
	static vector<vector<float> > gradients;
	static vector<vector<bool> > colorMasks;
	static vector<vector<int> > neighborRelations;
	static vector<int> labelIndices;
	
	TextureManager() {}
	
	static float computePlaneSmoothnessForMRF(int pix1, int pix2, int label1, int label2);
	vector<int> detectAndFillOcclusions(int noImages, int noPixels, int refIndex, int width, int height, vector<int> &pixels, vector<vector<Color> > &colors_, vector<vector<bool> > &colorMasks_, 
								 vector<int> &initialLabeling, vector<Color> &boundaryPixels, bool smoothBoundary);
};

#endif