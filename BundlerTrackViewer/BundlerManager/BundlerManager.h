/*
 *  BundlerManager.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/15/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _BUNDLER_MANAGER_H
#define _BUNDLER_MANAGER_H

#include <string>
#include "LineDetector.h"
#include "PerspectiveCamera.h"
#include "KeypointExtractor.h"
#include "KeypointMatcher.h"
#include "GeometryManager.h"
#include "BundleAdjustment.h"
#include "ImageRectification.h"
#include "RepetitionFinder.h"
#include "TranslationalGrid3D.h"
#include "PointCloud.h"
#include "GraphWrapper.h"
#include "Common.h"

class BundlerOutputParser;

using namespace std;

struct imagePatchVisibility
{
	int image;
	int visiblePatchNo;
	
	imagePatchVisibility() 
	{
		image = -1;
		visiblePatchNo = 0;
	}
	
	bool operator<(const imagePatchVisibility &rhs) const
	{
		if(visiblePatchNo > rhs.visiblePatchNo)
			return true;
		else
			return false;
	}
};

struct scaleEstimate
{
	int imgIndex1;
	int imgIndex2;
	int planeIndex1;
	int planeIndex2;
	float scale;
	int noInliers;
	
	bool operator<(const scaleEstimate &rhs) const
	{
		if(noInliers > rhs.noInliers)
			return true;
		else
			return false;
	}
};

struct gridShift
{
	int imageIndex1;
	int imageIndex2;
	int templateId;
	int planeIndex1;
	int planeIndex2;
	int gridIndex1;
	int gridIndex2;
	int shiftCol;
	int shiftRow;
	
	gridShift()
	{
		shiftCol = 0;
		shiftRow = 0;
	}
	
	bool operator==(const gridShift &rhs) const
	{
		if(rhs.imageIndex1 == imageIndex1 && rhs.imageIndex2 == imageIndex2 &&
			rhs.planeIndex1 == planeIndex1 && rhs.planeIndex2 == planeIndex2 &&
			rhs.gridIndex1 == gridIndex1 && rhs.gridIndex2 == gridIndex2 &&
			rhs.shiftCol == shiftCol && rhs.shiftRow == shiftRow)
				return true;
		else
			return false;
	}
};

struct imageTransformation
{
	int ind1;
	int ind2;
	bool gridMatch;
	vector<gridShift> gridShifts;
	float scale;
	Vec2f translation;
	Matrix3f fundMatrix;
	float score;
	int group;
	bool used;
	bool addedToGraph;
	
	imageTransformation()
	{
		ind1 = -1;
		ind2 = -1;
		gridMatch = false;
		scale = 1.0;
		translation = Vec2f(0.0, 0.0);
		fundMatrix.setIdentity();
		group = -1;
		used = false;
		addedToGraph = false;
		gridShifts.clear();
	}
	
	bool operator<(const imageTransformation &rhs) const
	{
		if(score > rhs.score)
			return true;
		else
			return false;
	}
	
	bool operator==(const imageTransformation &rhs) const
	{
		int g = 0;
		for(; g<rhs.gridShifts.size(); g++)
		{
			if(find(gridShifts.begin(), gridShifts.end(),rhs.gridShifts[g]) == gridShifts.end())
				break;
		}
		if(g == rhs.gridShifts.size())
			return true;
		else 
			return false;
	}
};

struct registeredImages {
	int templateId;
	vector<int> imageIndices;
	vector<int> planeIndices;
	vector<int> gridIndices;
	vector<Vec2i> shifts;
	repetitionGrid currentGrid;
	float score;
};

class BundlerManager
{
public:
	BundlerManager(string _mainDirectory);
	
	void setMainDirectory(string _mainDirectory) {mainDirectory = _mainDirectory; }
	void addNewImagePath(string path);
	void addNewMaskImagePath(string path) {maskFileNames.push_back(path);}
	
	Img* RectifyImage(int index, int& noPlanes);
	Img* getRectifiedImage(int imgIndex, int planeIndex);
	void computeEdges();
	
	void estimateScales();
	void findScaleOrder(int index, int planeIndex, vector<float> &scales);
	void readAllGridMatches();
	void readGridMatches(int index);
	void findMatches(int index, int rectificationInd, Vec2f topLeft, int w, int h);
	//void findPairwiseAlignments(GraphWrapper &gw, float graphEdgeThreshold, vector<vector<vector<imageTransformation > > > &imageAlignments, 
	//							vector<vector<vector<int > > > &inliers);
	//void findInconsistentCycles(vector<Cycle> &cycles, vector<vector<vector<imageTransformation > > > &imageAlignments, vector<vector<vector<int > > > &inliers);
	
	//int registerRepeatingRegions();
	//float evaluateImgToComponentAlignment(int refIndex, int neighIndex, int gridIndex, vector<int> component, vector<Vec2i> componentShifts, int shiftCol, int shiftRow);
	//int evaluateImgToComponentAlignment(int refIndex, int neighIndex, int gridIndex, vector<int> refComponent, vector<int> component, 
	//									vector<Vec2i> refComponentShifts, vector<Vec2i> componentShifts, int shiftCol, int shiftRow);
	//void findCompatibleMatches(vector<keyPairMatch> &tentativeMatches, 
	//						   vector<keyPair> &keyPairs1, vector<keyPair> &keyPairs2,
	//						   vector<keypt_t> &keyInfo1, vector<keypt_t> &keyInfo2,
	//						   vector<KeypointMatch> &bestMatches,
	//						   float closenessThresh, float distortionThresh);
	std::vector<std::vector<Vec2f> > getRepetitions(int index); 
	std::vector<std::vector<Vec2f> > getRefinedRepetitions(int index);
	std::vector<Color> getRepetitionColors(int index);
	void createRepetitionMask(int index);
	
	//std::vector<std::vector<Vec2f> > getGridStrokes(int index);
	//void projectGridSketchToImages(int imgIndex, vector<Vec2f> &endpoints);
	//void projectOccluderMaskToImages(int imageIndex);
	//void projectOccluderSketchToImages(int imageIndex, vector<Vec2f> &endpoints, int propagationType);
	//void removeOccluderUsingRepetitions(int imageIndex, vector<Vec2f> &endpoints);
	//void projectTextureSketchToImages(int imageIndex, vector<Vec2f> &endpoints);
	//void assignDepthToOcclusionArea(int minX, int minY, int maxX, int maxY, int planeIndex, int imageIndex, vector<int> &pixelIndices, vector<Vec3f> &positions);
	//void propagateOcclusionUsingDenseReconstruction(int minX, int maxX, int minY, int maxY, int imageIndex, int planeIndex, vector<int> &pixelIndices,  vector<Vec3f> &positions);
	//void warpImageToReferenceImage(int refIndex, int imageIndex, int planeIndex);
	//void fillOcclusionArea(int imageIndex, vector<Vec2f> &endpoints, vector<int> &pixelIndices, int noPixels);
	
	void ComputeAllFeaturePoints();
	void ReadFeaturePoints(int imageIndex);
	void ReadFeaturePointsWithRepetition(int imageIndex);
	void ReadAllFeaturePoints();
	void ReadAllFeaturePointsWithRepetition();
	void findUniqueFeaturePoints();
	void getFeaturePointInfo(int imageIndex, keypt_t* &featurePointInfo, int &noFeatures);
	void addFeaturePoint(int imageIndex, Vec2f pt);
	void saveFeaturePoints();
	
	void matchFeaturePoints(bool useMask);
	void buildManualCorrespondences();
	int getNoMatches(int imageIndex1, int imageIndex2);
	void getMatchingPoints(int imageIndex1, int imageIndex2, int matchIndex, Vec2f &pt1, Vec2f &pt2, Vec2f &epiStart, Vec2f &epiEnd, bool &valid);
	void getAllMatchingPoints(int imageIndex1, int imageIndex2, vector<Vec2f> &pt1, vector<Vec2f> &pt2, vector<bool> &valid);
	void setMatchValidity(int imageIndex1, int imageIndex2, int matchIndex, bool validity);
	void blendImagesWithMatches(int imageInd1, int imageInd2);
	
	void readMatches();
	void writeMatches();
	void saveManualFeatures();
	
	void computeTracks();
	int getNoOfTrack();
	int getNoOfCorrespondencesInTrack(int trackIndex);
	Correspondence getCorrespondenceInTrack(int trackIndex, int imageIndex);
	void formMatchesFromTracks();
	
	void runBundleAdjustment();
	void writeCalibrationResults();
	void readCalibrationResults();
	void generateCameraFrustums();
	void read3DGrids();
	
	float* getProjectedPoints(int &noPoints, int imageIndex);
	Color* getProjectedPointColors(int &noPoints, int imageIndex);
	
	//void cleanDensePointCloud();
	
	//void changeGridCount(int imageIndex, int gridIndex, int rowCount, int colCount);
	
	void bringToSameCoordinateFrame(vector<Vec2f> imageCoordinates1, vector<Vec2f> imageCoordinates2, int cam1, int cam2);
	
	void registerImages();
	void convertAlignmentsToMatches();
	void findCycles(GraphWrapper &gw,vector<vector<vector<imageTransformation> > > &imageAlignments);
	void findCandidateAlignments(GraphWrapper &gw,vector<vector<vector<imageTransformation> > > &imageAlignments);
	void evaluateAllPossibleAlignments(int refIndex, int neighIndex, int templateId, int refPlaneIndex, int neighPlaneIndex, int refGridIndex, int neighGridIndex, 
									   int &col, int &row, int &matches, vector<KeypointMatch> &tmpMatches);
	float evaluateImgToImgGridAlignment(int refIndex, int neighIndex, int templateId, int refPlaneIndex, int neighPlaneIndex, int refGridIndex, int neighGridIndex, 
										int shiftCol, int shiftRow, vector<KeypointMatch> &tmpMatches);
	float findAlignmentByDiscardingGrids(int refIndex, int neighIndex, int templateId, int refPlaneIndex, int neighPlaneIndex, int refGridIndex, int neighGridIndex, 
										 int &shiftCol, int &shiftRow, vector<KeypointMatch> &tmpMatches);
	float findMatchesForAlignment(int refIndex, int neighIndex, int templateId, int refPlaneIndex, int neighPlaneIndex, int refGridIndex, int neighGridIndex, 
								  int shiftCol, int shiftRow, vector<KeypointMatch> &tmpMatches);
	
	void fitPlanesTo3DData();
	void rectifyWith3DPlane();
	void updateGridCells();
	bool form3DGrids();
	
private:
	string mainDirectory;
	string bundlerOutputPath;
	vector<string> imageFileNames;
	vector<Vec2i> imageSizes;
	vector<Matrix3f> intrinsicMatrices;
	vector<string> maskFileNames;
	vector<PerspectiveCamera> cameras;
	
	//rectification
	vector<int> noRectifiedImages;
	vector<vector<Img*> > rectifiedImages;
	vector<vector<Matrix3f> > rectifyingHomographies;
	vector<Vec2f> vanishingPointDist; 
	ImageRectification *imgRect;
	
	//edges
	vector<vector<Line> > lines;
	LineDetector *lineDetector;
	
	//keypoints
	KeypointExtractor *kpExtractor;
	vector<vector<short> > keys;
	vector<vector<keypt_t> > keyInfo;
	vector<int> noFeaturePoints;
	bool keySpaceAllocated;
	vector<vector<keyPair> > keyPairs;
	
	//matches
	KeypointMatcher *kpMatcher;
	vector<vector<vector<KeypointMatch> > > matches;
	vector<vector<Matrix3f> > fundMatrices;
	
	//tracks
	std::vector<CorrespondenceTrack> tracks;
	
	//repetitions
	vector<vector<scaleEstimate> > scaleEstimates;
	vector<vector<int> > gridMask;
	vector<vector<int> > shiftCols;
	vector<vector<int> > shiftRows;
	repetitionGrid finalGrid;
	vector<vector<vector<Vec2i> > > elementSizes;//for each template->for each image->for each plane
	vector<Color> repetitionMaskColors;
	vector<vector<vector<vector<vector<vector<Vec2f> > > > > > repetitions;//for each template->for each image->for each plane->for each grid
	vector<vector<vector<vector<vector<Vec2f> > > > > refinedRepetitions;//for each template->for each image->for each plane->for each grid
	vector<vector<vector<vector<repetitionGrid> > > > repeatingGrids; //for each template->for each image->for each plane->each grid
	vector<TranslationalGrid3D> grids3D;
	vector<imageTransformation> usedTransformations;
	vector<vector<vector<pair<int, int> > > > imageGridIndices;
	
	//bundle adjustment
	std::vector<std::vector<Transformation> > transformations;
	GeometryManager *geoManager;
	BundlerOutputParser *boParser;
	vector<Matrix4f> camMatrices;
	vector<Vec3f> points3D;
	vector<vector<Vec2f> > pointProjections;
	vector<vector<Color> > pointProjectionColors;
	
	//edits
	vector<vector<vector<Vec2f> > > gridStrokes;
	
	//dense and sparse reconstruction (if exists)
	PointCloud *densePC;
	PointCloud *sparsePC;
	
	vector<Plane> planes;
	vector<vector<Vec3f> > planeInliers; 
};

#endif
