/*
 *  MRFDepthLabeling.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 4/16/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include <algorithm>
#include <sstream>
#include "MRFDepthLabeling.h"
#include "EpipolarGeometry.h"
#include "../MRFUtils/MRFController.h"

#define COARSE_LEVEL_WEIGHT 0.5
#define GRAD_THRESH 0.03
#define GRAD_SCALE 2.0

vector<float> MRFDepthLabeling::vCues = vector<float> ();
vector<float> MRFDepthLabeling::hCues = vector<float> ();
vector<int> MRFDepthLabeling::labelIndices = vector<int> ();
vector<vector<int> > MRFDepthLabeling::neighborRelations = vector<vector<int> >();
int MRFDepthLabeling::smoothScale = 1.0;
int MRFDepthLabeling::currentWidth = 0;
int MRFDepthLabeling::currentHeight = 0;

MRFDepthLabeling::MRFDepthLabeling(vector<string> &imageFilenames_, vector<PerspectiveCamera> &cameras_)
{
	for(int i=0; i<imageFilenames_.size(); i++)
		imageFilenames.push_back(imageFilenames_[i]);
	for(int i=0; i<cameras_.size(); i++)
		cameras.push_back(cameras_[i]);
}

float MRFDepthLabeling::computeSmoothnessForMRF(int pix1, int pix2, int label1, int label2)
{
	double lambda     = 0.25*smoothScale;
    double smoothExp = 2;
    double smoothMax = 5.0;
	
	if(find(neighborRelations[pix1].begin(), neighborRelations[pix1].end(), pix2) == neighborRelations[pix1].end())
	{
		printf("Error in neighboring between variable %d and %d!\n", pix1, pix2);
		return 0.0;
	}

	int pixelIndex1 = MRFDepthLabeling::labelIndices[pix1];
	int pixelIndex2 = MRFDepthLabeling::labelIndices[pix2];
	
	int x1 = pixelIndex1 % currentWidth;
	int y1 = pixelIndex1 / currentWidth;
	
	int x2 = pixelIndex2 % currentWidth;
	int y2 = pixelIndex2 / currentWidth;
	
	if(x1 == x2)
	{
		//vertical
		if(y1<y2)
		{
			return lambda * std::min(smoothMax, pow(abs(label1-label2), smoothExp))*vCues[pixelIndex1];
		}
		else
		{
			return lambda * std::min(smoothMax, pow(abs(label1-label2), smoothExp))*vCues[pixelIndex2];
		}
	}
	else if(y1 == y2)
	{
		//horizontal
		if(x1 < x2)
		{
			return lambda * std::min(smoothMax, pow(abs(label1-label2), smoothExp))*hCues[pixelIndex1];
		}
		else
		{
			return lambda * std::min(smoothMax, pow(abs(label1-label2), smoothExp))*hCues[pixelIndex2];
		}
	}
	else
	{
		printf("Error in MRF optimization, non-grid neighbors!\n");
	}
}

void MRFDepthLabeling::findNeighboringViews(int index, vector<NeighborView>& neighboringViews) 
{
	// Find images with viewing angle constraints 
	//find median dist between cameras
	float medianDistBetweenCam = 0;
	int count=0;
	for(int i=0; i<cameras.size(); i++)
	{
		for(int j=i+1;j<cameras.size();j++)
		{
			medianDistBetweenCam += (cameras[i].getCenter() - cameras[j].getCenter()).length();
			count++;
		}
	}
	medianDistBetweenCam /= float(count);
	
	
	
	float viewAngleThresh = 15;
	neighboringViews.clear();
	
	float angle;
	int noCandidates;
	
	Vec4f ray = cameras[index].getOAxis();
	ray[3] = 0.0f;
	
	noCandidates = cameras.size();
	for(int i=0; i<noCandidates; i++)
	{
		if(i==index)
			continue;
		
		Vec4f rayNeighbor = cameras[i].getOAxis();
		rayNeighbor[3] = 0.0f;
		
		angle = dot(ray, rayNeighbor);
		
		if(angle < cos(viewAngleThresh*M_PI/180.0))
			continue;
		
		float dist = (cameras[index].getCenter() - cameras[i].getCenter()).length();
		
		if(dist > 2*medianDistBetweenCam || dist < 0.005*medianDistBetweenCam)
			continue;
		
		NeighborView nv(i, angle);
		neighboringViews.push_back(nv);
	}
	sort(neighboringViews.begin(), neighboringViews.end());
}

int MRFDepthLabeling::resizePixelIndices(vector<int> &oldPixelIndices, vector<int> &newPixelIndices, int widthPrev, int heightPrev, int width, int height)
{
	float ratio_x = float(widthPrev) / float(width);
	float ratio_y = float(heightPrev) / float(height);
	newPixelIndices.resize(width*height, -1);
	
	MRFDepthLabeling::labelIndices.clear();
	
	int noLabels = 0;
	for (int y = 0; y < height; y++) 
	{
		for (int x = 0; x < width; x++) 
		{
			float xf, yf;
			int xi, yi;
			
			xf = float(x) * ratio_x;
			yf = float(y) * ratio_y;
			xf = min(max(xf, 0.0f), float(widthPrev-1));
			yf = min(max(yf, 0.0f), float(heightPrev-1));
			
			xi = floor(xf);
			yi = floor(yf);
			if(oldPixelIndices[yi*(widthPrev)+xi] != -1)
			{
				labelIndices.push_back(y*width+x);
				newPixelIndices[y*width+x] = noLabels;
				noLabels+=1;
			}
		}
	}
	return noLabels;
}

void MRFDepthLabeling::resizeLabels(vector<int> &prevVector, vector<bool> &prevVectorFlag, int widthPrev, int heightPrev, int width, int height, int resPrev, int res)
{
	if (prevVector.size() > 0) 
	{
    	float ratio_x = float(widthPrev) / float(width);
    	float ratio_y = float(heightPrev) / float(height);
    	float ratio_r = float(resPrev) / float(res);
		vector<float> newVector;
		vector<bool> newVectorFlag;
		newVector.resize(width*height, 0.0);
		newVectorFlag.resize(width*height, false);
		
#pragma omp parallel for
		for (int y = 0; y < height; y++) 
		{
    	    for (int x = 0; x < width; x++) 
			{
				float xf, yf;
				int xi, yi;
				
				xf = float(x) * ratio_x;
				yf = float(y) * ratio_y;
				xf = min(max(xf, 0.0f), float(widthPrev-1));
				yf = min(max(yf, 0.0f), float(heightPrev-1));
				
				xi = floor(xf);
				yi = floor(yf);
				newVectorFlag[y*width+x] = prevVectorFlag[yi*widthPrev+xi];
				newVector[y*width+x] = prevVector[yi*widthPrev+xi];
				
			}
    	}
		prevVector.clear();
		prevVector.assign(newVector.begin(), newVector.end());
		prevVectorFlag.clear();
		prevVectorFlag.assign(newVectorFlag.begin(), newVectorFlag.end());
    } 
	else 
	{
		prevVector.resize(width*height, 0.0);
		prevVectorFlag.resize(width*height, false);
    }
}


void MRFDepthLabeling::resizeEnergy(vector<float> &prevVector, int widthPrev, int heightPrev, int width, int height, int resPrev, int res)
{
	if(prevVector.size() > 0)
	{
		float ratio_x = float(widthPrev) / float(width);
    	float ratio_y = float(heightPrev) / float(height);
    	float ratio_r = float(resPrev) / float(res);
		vector<float> newVector;
		newVector.resize(width*height*res, 0.0);
		for (int y = 0; y < height; y++) 
		{
    	    for (int x = 0; x < width; x++) 
			{
				for (int r = 0; r < res; r++) 
				{
					
					float xf, yf, rf;
					int xi, yi, ri;
					
					xf = float(x) * ratio_x;
					yf = float(y) * ratio_y;
					rf = float(r) * ratio_r;
					
					xf = min(max(xf, 0.0f), float(widthPrev-1));
					yf = min(max(yf, 0.0f), float(heightPrev-1));
					rf = min(max(rf, 0.0f), float(resPrev-1));
					
					xi = min(int(xf), widthPrev-2);
					yi = min(int(yf), heightPrev-2);
					ri = min(int(rf), resPrev-2);
					
					float val1 = (1.0-xf+xi)*((1.0-rf+ri)*prevVector[(yi*widthPrev+xi)*resPrev+ri] + (rf-ri)*prevVector[(yi*widthPrev+xi)*resPrev+ri+1]) + 
					(xf-xi)*((1.0-rf+ri)*prevVector[(yi*widthPrev+xi+1)*resPrev+ri] + (rf-ri)*prevVector[(yi*widthPrev+xi+1)*resPrev+ri+1]);
					float val2 = (1.0-xf+xi)*((1.0-rf+ri)*prevVector[((yi+1)*widthPrev+xi)*resPrev+ri] + (rf-ri)*prevVector[((yi+1)*widthPrev+xi)*resPrev+ri+1]) + 
					(xf-xi)*((1.0-rf+ri)*prevVector[((yi+1)*widthPrev+xi+1)*resPrev+ri] + (rf-ri)*prevVector[((yi+1)*widthPrev+xi+1)*resPrev+ri+1]);
					
					newVector[(y*width+x)*res + r] = (1.0-yf+yi)*val1 + (yf-yi)*val2;
					
				}
    	    }
    	}
		prevVector.clear();
		prevVector.assign(newVector.begin(), newVector.end());
	}
	else 
	{
		prevVector.resize(width*height*res, -1.0);
	}
	
}

void MRFDepthLabeling::computeGradientSensitiveMrfCues(Img &refImg, vector<int> &currentPixelIndices, int minX, int maxX, int minY, int maxY, vector<float> &vCues, vector<float> &hCues)
{
    int width = maxX - minX;
	int height = maxY - minY;
	
    hCues.resize(width*height, 0.0); 
	vCues.resize(width*height, 0.0);
	
	for(int x=minX;x<maxX;x++)
	{
		for(int y=minY;y<maxY;y++)
		{
			if(currentPixelIndices[(y-minY)*width+(x-minX)] == -1)
			{
				continue;
			}
			
			if(currentPixelIndices[(y-minY)*width+((x+(x<maxX-1))-minX)] != -1)
			{
				Color current = refImg(x,y);
				Color nextHor = refImg(x+(x<maxX-1), y);
				float sx = 0.0;
				for(int i=0; i<3; i++)
				{
					sx += (current[i]-nextHor[i])*(current[i]-nextHor[i]);
				}
				if(sx < GRAD_THRESH)
					hCues[(y-minY)*width+(x-minX)] = GRAD_SCALE;
				else
					hCues[(y-minY)*width+(x-minX)] = 1.0;
			}
			if(currentPixelIndices[((y+(y<maxY-1))-minY)*width+(x-minX)] != -1)
			{
				Color current = refImg(x,y);
				Color nextVer = refImg(x, y+(y<maxY-1));
				float sy = 0.0;
				for(int i=0; i<3; i++)
				{
					sy += (current[i]-nextVer[i])*(current[i]-nextVer[i]);
				}
				if(sy < GRAD_THRESH)
					vCues[(y-minY)*width+(x-minX)] = GRAD_SCALE;
				else
					vCues[(y-minY)*width+(x-minX)] = 1.0;
			}
		}
	}
}

void MRFDepthLabeling::computeDepthMapForImagePair(int indexRef, int indexNeighbor, int level, Img& refImg, Img& neighborImg, int minX, int maxX, int minY, int maxY,
												   vector<int> &pixelIndices, vector<float> &nccEnergies, int levelRes, bool limitSearch,
												   float depthOffset, vector<int> &prevDepths, vector<bool> &prevDepthAssigned)
{
	int width = (maxX-minX);
	int height = (maxY-minY);
	
	vector<float> nccScores;
	nccScores.resize(levelRes*width*height);
	
	clock_t begin, end;
	
	EpipolarGeometry eg;
	
	begin = clock();
	
#pragma omp parallel for
	for(int y=minY; y<maxY; y++)
	{
		for(int x=minX; x<maxX; x++)
		{
			if(pixelIndices[(y-minY)*width+(x-minX)] == -1)
				continue;
			
			Vec3f pixelPos(x, y, 1.0);
			float res;
			int minLabel = 0;
			float minD, maxD;
			
			//assign min and max depths
			if(prevDepths.size() > 0)
			{
				int pixelIndex = (y-minY)*width+(x-minX);				
				if(prevDepths[pixelIndex] != 0 && limitSearch)
				{
					minD = (minDepth + ((maxDepth - minDepth)/levelRes * prevDepths[pixelIndex])) - depthOffset;
					maxD = (minDepth + ((maxDepth - minDepth)/levelRes * prevDepths[pixelIndex])) + depthOffset;
						
					if(minD < minDepth)
						minD = minDepth;
					if(maxD > maxDepth)
						maxD = maxDepth;
						
					//how many depth intervals do we need to use between min and max depth
					if(maxD-minD < 2*depthOffset)
						res = (depthRes * (maxD-minD)) / (2.0*depthOffset);
					else
						res = depthRes;
				}			
				
				else
				{
					minD = minDepth;
					maxD = maxDepth;
					res = levelRes;
				}
				
				minLabel = ((minD - minDepth)*levelRes) / (maxDepth-minDepth);
				
				eg.findSamplesOnEpipolarLine(pixelPos, pixelIndex, cameras[indexRef], cameras[indexNeighbor], refImg, neighborImg, level,
											 minDepth, maxDepth, res, nccEnergies, minLabel, levelRes);
			}
			
		}
	}
	end = clock();
	printf("NCC at level %d:%f\n", level, double((end-begin)*1000.0/CLOCKS_PER_SEC));
}

void MRFDepthLabeling::computeDepthMap(int minX, int maxX, int minY, int maxY, vector<int> &pixelIndices, float minDepth_, float maxDepth_, int refIndex, vector<float> &depths)
{
	vector<NeighborView> neighboringViews;
	findNeighboringViews(refIndex, neighboringViews);
	//neighboringViews[0].index = 10;
	//neighboringViews[1].index = 11;
	//neighboringViews[2].index = 13;
	int noNeigborViews = neighboringViews.size();
	
	//int maxNoNeighborViews = 3;
	//if(noNeigborViews > maxNoNeighborViews)
	//	noNeigborViews = maxNoNeighborViews;
	
	int startLevel = 2;
	int endLevel = 0;
	bool limitSearch = false;
	
	minDepth = minDepth_;
	maxDepth = maxDepth_;
	
	PerspectiveCamera refCam = cameras[refIndex];
	Vec3f camCenter = shrink4To3(refCam.getCenter());
	
	std::vector<int> prevDepths;
	std::vector<bool> prevDepthAssigned;
	prevDepths.clear(); prevDepthAssigned.clear();
	
	//mrf
	std::vector<float> nccEnergies;
	std::vector<float> dataEnergies;
	std::vector<int> finalDepths;
	nccEnergies.clear();
	dataEnergies.clear();
	
	int noPixels;
	int widthPrev, heightPrev;
	float depthOffset = (maxDepth - minDepth) / 2.0;
	int levelRes = 32;
	int depthRes = 32;
	bool coarseLevel = true;
	
	string filename("/Users/ceylan/Desktop/depthmap");
	widthPrev = 0; heightPrev = 0;	
	
	clock_t begin, end;
	
	for(int j=startLevel; j>=endLevel; j--)
	{
		begin = clock();
		
		if(j<startLevel)
			coarseLevel = false;
		
		int currentMinX = minX / pow(2.0, j);
		int currentMaxX = maxX / pow(2.0, j);
		int currentMinY = minY / pow(2.0, j);
		int currentMaxY = maxY / pow(2.0, j);
		
		currentWidth = currentMaxX - currentMinX;
		currentHeight = currentMaxY - currentMinY;
		noPixels = currentWidth*currentHeight;
		
		if(j<startLevel)
		{
			widthPrev = (maxX-minX) / pow(2.0, j+1);
			heightPrev = (maxY-minY) / pow(2.0, j+1);
		}
		
		vector<int> currentPixelIndices;
		
		int noVar  = resizePixelIndices(pixelIndices, currentPixelIndices, (maxX-minX), (maxY-minY), currentWidth, currentHeight);
		
		resizeLabels(prevDepths, prevDepthAssigned, widthPrev, heightPrev, currentWidth, currentHeight, levelRes, depthRes*pow(2.0, startLevel-j));
		resizeEnergy(dataEnergies, widthPrev, heightPrev, currentWidth, currentHeight, levelRes, depthRes*pow(2.0, startLevel-j));
		
		//fill neighboring relations
		neighborRelations.clear();
		for(int y=0; y<currentHeight; y++)
		{
			for(int x=0; x<currentWidth; x++)
			{
				if(currentPixelIndices[y*currentWidth+x] != -1)
				{
					vector<int> neighbors;
					
					//printf("neighbors for variable:%d:", currentPixelIndices[y*currentWidth+x]);
					if(x-1>0)
					{
						if(currentPixelIndices[y*currentWidth + (x-1)] != -1)
						{
							neighbors.push_back(currentPixelIndices[y*currentWidth + (x-1)]);
							//printf("%d ", currentPixelIndices[y*currentWidth + (x-1)]);
						}
					}
					if(x+1<currentWidth)
					{
						if(currentPixelIndices[y*currentWidth + (x+1)] != -1)
						{
							neighbors.push_back(currentPixelIndices[y*currentWidth + (x+1)]);
							//printf("%d ", currentPixelIndices[y*currentWidth + (x+1)]);
						}
					}
					if(y-1>0)
					{
						if(currentPixelIndices[(y-1)*currentWidth + x] != -1)
						{
							neighbors.push_back(currentPixelIndices[(y-1)*currentWidth + x]);
							//printf("%d ", currentPixelIndices[(y-1)*currentWidth + x]);
						}
					}
					if(y+1<currentHeight)
					{
						if(currentPixelIndices[(y+1)*currentWidth + x] != -1)
						{
							neighbors.push_back(currentPixelIndices[(y+1)*currentWidth + x]);
							//printf("%d ", currentPixelIndices[(y+1)*currentWidth + x]);
						}
					}
					//printf("\n");
					MRFDepthLabeling::neighborRelations.push_back(neighbors);
				}
				
			}
		}

		if(limitSearch)
			depthOffset = ((maxDepth - minDepth) / 2.0) / pow(2.0, startLevel - j); // depth offset to constrain the search
		levelRes = depthRes*pow(2.0, startLevel-j); // total number of labels for this level
		nccEnergies.clear();
		nccEnergies.resize(noPixels*levelRes, 0.0);
		
		Img refImg;
		printf("filename:%s\n", imageFilenames[refIndex].c_str());
		refImg.read(imageFilenames[refIndex]);
		
		//resize the image
		for(int i=0; i<j; i++)
		{
			refImg.downScaleAndGaussianSmoothImage();
		}
		
		for(int i=0; i<noNeigborViews; i++)
		{
			Img neighborImg;
			printf("filename:%s\n", imageFilenames[neighboringViews[i].index].c_str());
			neighborImg.read(imageFilenames[neighboringViews[i].index]);
			
			for(int level=0; level<j; level++)
			{
				neighborImg.downScaleAndGaussianSmoothImage();
			}
			
			computeDepthMapForImagePair(refIndex, neighboringViews[i].index, j, refImg, neighborImg, currentMinX, currentMaxX, currentMinY, currentMaxY, 
										currentPixelIndices, nccEnergies, levelRes, limitSearch, depthOffset,  prevDepths, prevDepthAssigned);
			
		}
		
		begin = clock();
		
		vector<float> tmpDataEnergies;
		tmpDataEnergies.resize(noVar*levelRes);
		
		vector<int> initialLabels;
		initialLabels.resize(noVar);
		
#pragma omp parallel for
		for(int i=0; i<dataEnergies.size(); i++)
		{
			int pix = i / levelRes;
			int label = i % levelRes;
			
			if(j<startLevel)
			{
				dataEnergies[i] = COARSE_LEVEL_WEIGHT*dataEnergies[i] + (1.0-COARSE_LEVEL_WEIGHT)*(1.0-nccEnergies[i]);
			}
			else
			{
				dataEnergies[i] = 1.0 - nccEnergies[i];
			}
			if(currentPixelIndices[pix] != -1)
			{
				tmpDataEnergies[currentPixelIndices[pix]*levelRes + label] = dataEnergies[i];
				initialLabels[currentPixelIndices[pix]] = prevDepths[pix];
			}
		}
		
		computeGradientSensitiveMrfCues(refImg, currentPixelIndices, currentMinX, currentMaxX, currentMinY, currentMaxY, MRFDepthLabeling::vCues, MRFDepthLabeling::hCues);
		
		MRFController *mrfController = new MRFController();
		
		//mrfController->performMRFWithGraphCutExpansion(tmpDataEnergies, MRFDepthLabeling::computeSmoothnessForMRF, prevDepths, finalDepths, currentWidth, currentHeight, levelRes);
	
		mrfController->performMRFWithGraphCutExpansion(tmpDataEnergies, MRFDepthLabeling::neighborRelations, MRFDepthLabeling::computeSmoothnessForMRF, 
													   initialLabels, finalDepths, noVar, levelRes);
		
		delete mrfController;
		
		prevDepths.clear();
		prevDepths.resize(noPixels);
		
		#pragma omp parallel for
		for(int x=0; x<currentWidth; x++)
		{
			for(int y=0; y<currentHeight; y++)
			{
				int labelIndex = currentPixelIndices[y*currentWidth+x];
				if(labelIndex == -1)
					prevDepths[y*currentWidth+x] = 0;
				else
					prevDepths[y*currentWidth+x] = finalDepths[labelIndex];
			}
		}
		
		//prevDepths.assign(finalDepths.begin(), finalDepths.end());
		
		float e = 0.0;
		vector<float> colors;
		colors.resize(noPixels, 0.0);
		
		for(int i=0;i<noPixels;i++)
		{
			if(currentPixelIndices[i] != -1)
			{
				colors[i] = 1.0 - (float(finalDepths[currentPixelIndices[i]]) / float(levelRes));
			}
		}
		
		if(j==endLevel)
		{
			depths.resize(noPixels, 0);
			float delta = (maxDepth - minDepth) / float(levelRes);
			for(int i=0;i<noPixels;i++)
			{
				if(currentPixelIndices[i] != -1)
				{
					depths[i] = minDepth +  finalDepths[currentPixelIndices[i]]*delta;
				}
			}
		}
		
		end = clock();
			
		printf("MRF at level %d:%f\n", j, double((end-begin)*1000.0/CLOCKS_PER_SEC));
		
		//save result
		Img *depthMap = new Img(currentWidth, currentHeight, colors, true);
		stringstream ss;
		ss << filename << j << ".png";
		string s;
		ss >> s;
		depthMap->write(s);	
	}
	
	if(endLevel != 0)
	{
		resizeLabels(prevDepths, prevDepthAssigned, currentWidth, currentHeight, maxX-minX, maxY-minY, levelRes, depthRes*pow(2.0, startLevel-endLevel));
		noPixels = (maxX-minX)*(maxY-minY);
		vector<float> colors;
		colors.resize(noPixels, 0.0);
		depths.resize(noPixels, 0.0);
		float delta = (maxDepth - minDepth) / float(levelRes);
		
		for(int i=0;i<noPixels;i++)
		{
			depths[i] = minDepth +  prevDepths[i]*delta;
			colors[i] = 1.0 - ((depths[i]-minDepth)/(maxDepth-minDepth));
		}
		
		Img *depthMap = new Img( maxX-minX, maxY-minY, colors, true);
		stringstream ss;
		ss << filename << 0 << ".png";
		string s;
		ss >> s;
		depthMap->write(s);	
	}
}