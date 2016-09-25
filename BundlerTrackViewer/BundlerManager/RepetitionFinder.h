/*
 *  RepetitionFinder.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 12/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _REPETITION_FINDER_H
#define _REPETITION_FINDER_H

#include "Photo.h"
#include "Image.h"
#include "Common.h"

#include <vector>

using namespace std;

struct gridTransformation 
{
	int startIndex;
	int endIndex;
	int transformationId;
	Vec2f transformation;
};

struct repetitionCell 
{
	Vec2f center;
	Vec2i size;
	int rowNumber;
	int colNumber;
	float descriptor;
	vector<float> pcaCoord;
	vector<vector<float> > allPcaCoord;
	int label;
	int visImages;
	bool temporary;
	
	repetitionCell(Vec2f c)
	{
		center = c;
		rowNumber = -1;
		colNumber = -1;
		descriptor = -1.0;
		visImages = 1;
		temporary = false;
	}
};

struct repetitionGrid
{
	vector<repetitionCell> gridCells;
	Vec2f xTrans, yTrans;
	int numColumns;
	int numRows;
};

struct gridRegistrationCandidate
{
	int colShift;
	int rowShift;
	float score;
	
	bool operator<(const gridRegistrationCandidate &rhs) const
	{
		if(score < rhs.score)
			return true;
		else
			return false;
	}
};

class RepetitionFinder
{
public:
	RepetitionFinder() {}
	
	Img getAvgTemplate(vector<Img *> &templateImg, vector<repetitionGrid> &matches, vector<Vec2i> &templateSizes, Vec2i orgTempSize);
	
	void performPCA(vector<Img *> &templateImg, vector<vector<repetitionGrid> > &matches, int gridIndex, vector<Vec2i> &templateSizes);
	
	float computeNCCScoreWithGrayScale(vector<Color> &refPatch, vector<Color> &neighborPatch);
	float computeNCCScoreWithGrayScale(Vec2f refPixel, Vec2f neighborPixel, Photo &refPhoto, Photo &neighPhoto, float windowSize, int level);
	int findBestMatch(Vec2f refPixel, vector<Vec2f> &candidateMatches, Photo &refPhoto, Photo &neighPhoto, int level);
	
	void templateMatchingWithSmallerFFTCorrelation(Img* templateImg, Img* sourceImg, vector<Vec2f> &matches, vector<float> &scores, float thresh);
	void templateMatchingWithFFTCorrelation(Img* templateImg, Img* sourceImg, vector<Vec2f> &matches, vector<float> &scores, float thresh);
	void templateMatchingWithFFTCorrelation(Img* templateImg, Img* sourceImg, Vec2i offset, Vec2i searchSize, vector<Vec2f> &matches, vector<float> &scores, float thresh);
	void findRepetitionInROI(Img *sourceImg, Img* templateImg, Vec2f upperCorner, Vec2f lowerCorner, vector<Vec2f> &matches, vector<float> &scores, float thresh);
	
	void findGridGeneratingVectors(vector<Vec2f> &matchingCenters, vector<float> &scores, float xThresh, float yThresh, 
								   vector<vector<Vec2f> > &finalGroups, vector<Vec2f> &finalTransformations);
	void segmentIntoRowsAndColumns(vector<Vec2f> &matchingCenters, vector<float> &scores, repetitionGrid &finalGrid, Vec2f gridVectorY, Vec2f gridVectorX, int w, int h);
	void fillMissingGridCells(repetitionGrid &grid, Img *templateImg, Img *sourceImg, float thresh);
	void fillMissingGridCellsTemporarily(repetitionGrid &grid, Vec2i size);
	
	int translateTransformationToGrid(repetitionGrid &refGrid, repetitionGrid &neighborGrid, float scale, Vec2f translation, 
									  float repThresh, int width, int height, int &colShift, int &rowShift);
	int translateTransformationToGrid(repetitionGrid &refGrid, repetitionGrid &neighborGrid, vector<Vec2f> &aPts, vector<Vec2f> &bPts, 
									   vector<int> &inliers, float scale, Vec2f translation, float repThresh, int width, int height, int &shiftCol, int &shiftRow);
	float translateTransformationToGridRANSAC(repetitionGrid &refGrid, repetitionGrid &neighborGrid, vector<Vec2f> &aPts, vector<Vec2f> &bPts, 
											  vector<float> &scale, vector<Vec2f> &translation, vector<vector<int> > &matches, 
											  float repThresh, int width, int height, vector<int> &shiftCol, vector<int> &shiftRow, vector<vector<int> > &inlierMatches,
											  int &bestCandidate);
	bool refineMatching(repetitionGrid &refGrid, repetitionGrid &neighborGrid, repetitionGrid &currentGrid, Vec2i neighborShift, 
						float repThresh, float scale, Vec2f &translate, float ratio, vector<int> &shiftCol, vector<int> &shiftRow);
	void findBestMatchingScore(repetitionGrid &refGrid, repetitionGrid &neighborGrid, float &score, float &ratio, int &shiftCol, int &shiftRow);
	repetitionGrid mergeGrids(repetitionGrid &grid1, repetitionGrid &grid2, int shiftCol, int shiftRow);
	
private:
	//cv::PCA *pca;
};

#endif