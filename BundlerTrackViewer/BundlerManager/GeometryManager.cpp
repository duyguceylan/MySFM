/*
 *  GeometryManager.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "GeometryManager.h"
#include "matrix/matrix.h"
#include <cmath>

GeometryManager::GeometryManager()
{
	eg = NULL;
	
	fMatrixMaxTrialNumber = 2048;
	fMatrixMinInlierNumber = 16;
	fMatrixThreshold = 9.0;
	
	homographyTrialNumber = 256;
	homographyThreshold = 6.0;
	homographyInlierNumber = 10;
}

bool GeometryManager::computeEpipolarGeometry(int imageIndex1, int imageIndex2, bool removeBadMatches,
											  vector<vector<keypt_t> > &keyInfo, vector<vector<vector<KeypointMatch> > > &matches, 
											  std::vector<std::vector<Transformation> > &transformations) 
{
    Matrix3f F;
    
	if(eg == NULL)
		eg = new EpipolarGeometry();
	
	std::vector<int> inliers;
	int numInliers;
	
	if(imageIndex1 < imageIndex2)
	{
		std::vector<KeypointMatch> &list = matches[imageIndex1][imageIndex2-imageIndex1-1];
		inliers  = eg->estimateFMatrix(keyInfo[imageIndex1], keyInfo[imageIndex2], list, fMatrixMaxTrialNumber, fMatrixThreshold, F);
		numInliers = (int) inliers.size();
		
		if (removeBadMatches) 
		{
			// Refine the matches
			std::vector<KeypointMatch> newMatchList;
			
			for (int i = 0; i < numInliers; i++) 
			{
				newMatchList.push_back(list[inliers[i]]);
			}
			
			list.clear();
			list = newMatchList;
		}
	}
	else
	{
		std::vector<KeypointMatch> &list = matches[imageIndex2][imageIndex1-imageIndex2-1];
		inliers =  eg->estimateFMatrix(keyInfo[imageIndex2], keyInfo[imageIndex1], list, fMatrixMaxTrialNumber, fMatrixThreshold, F);
		numInliers = (int) inliers.size();
		
		if (removeBadMatches) 
		{
			// Refine the matches
			std::vector<KeypointMatch> newMatchList;
			
			for (int i = 0; i < numInliers; i++) 
			{
				newMatchList.push_back(list[inliers[i]]);
			}
			
			list.clear();
			list = newMatchList;
		}
	}
	
	printf("numInliers:%d\n", numInliers);
    
    if (numInliers >= fMatrixMinInlierNumber) 
	{
		if(imageIndex1 < imageIndex2)
		{
			transformations[imageIndex2][imageIndex1].fMatrix = F.transpose();
			transformations[imageIndex1][imageIndex2].fMatrix = F;
		}
		else
		{
			transformations[imageIndex1][imageIndex2].fMatrix = F.transpose();
			transformations[imageIndex2][imageIndex1].fMatrix = F;
		}
		
		Matrix3f F_ = transformations[imageIndex1][imageIndex2].fMatrix;
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				printf("%f ", F_[j][i]);
			}
			printf("\n");
		}
		return true;
	} 
	else 
	{
		return false;
    }
}

void GeometryManager::computeEpipolarGeometry(std::vector<std::vector<Transformation> > &transformations, bool removeBadMatches, vector<vector<keypt_t> > &keyInfo, vector<vector<vector<KeypointMatch> > > &matches)
{
	unsigned int numImages = matches.size();
	
	for (unsigned int i = 0; i < numImages; i++) 
	{
        for(unsigned int iter=0; iter<matches[i].size(); iter++ )
		{
			bool connect = computeEpipolarGeometry(i, iter+i+1, removeBadMatches, keyInfo, matches, transformations);
			if(connect)
			{
				Matrix3f F_ = transformations[i][iter+i+1].fMatrix;
				for(int i=0; i<3; i++)
				{
					for(int j=0; j<3; j++)
					{
						printf("%f ", F_[j][i]);
					}
					printf("\n");
				}
			}
			else
			{
				//remove matches
				if(removeBadMatches)
				{
					matches[i][iter+i+1].clear();
				}
			}
		}
	}
	
}

/* Compute a transform between a given pair of images */
bool GeometryManager::computeTransform(int imageIndex1, int imageIndex2, bool removeBadMatches,
								  vector<vector<keypt_t> > &keyInfo, vector<vector<vector<KeypointMatch> > > &matches, 
								  std::vector<std::vector<Transformation> > &transformations)
{
	Matrix3f M;
	
	if(eg == NULL)
		eg = new EpipolarGeometry();
	
	std::vector<int> inliers;
	int numInliers;
	
	if(imageIndex1 < imageIndex2)
	{
		std::vector<KeypointMatch> &list = matches[imageIndex1][imageIndex2-imageIndex1-1];
		//inliers  = eg->estimateTransform(keyInfo[imageIndex1], keyInfo[imageIndex2], list, homographyTrialNumber, homographyThreshold, M);
		numInliers = (int) inliers.size();
		
		if (removeBadMatches) 
		{
			// Refine the matches
			std::vector<KeypointMatch> newMatchList;
			
			for (int i = 0; i < numInliers; i++) 
			{
				newMatchList.push_back(list[inliers[i]]);
			}
			
			list.clear();
			list = newMatchList;
		}
	}
	else
	{
		std::vector<KeypointMatch> &list = matches[imageIndex2][imageIndex1-imageIndex2-1];
		//inliers =  eg->estimateTransform(keyInfo[imageIndex2], keyInfo[imageIndex1], list, homographyTrialNumber, homographyThreshold, M);
		numInliers = (int) inliers.size();
		
		if (removeBadMatches) 
		{
			// Refine the matches
			std::vector<KeypointMatch> newMatchList;
			
			for (int i = 0; i < numInliers; i++) 
			{
				newMatchList.push_back(list[inliers[i]]);
			}
			
			list.clear();
			list = newMatchList;
		}
	}
	
	if (numInliers >= homographyInlierNumber) 
	{
		if(imageIndex1 < imageIndex2)
		{
			transformations[imageIndex2][imageIndex1].hMatrix = M.getInverseMatrix();
			transformations[imageIndex1][imageIndex2].hMatrix = M;
			
			transformations[imageIndex1][imageIndex2].inlierNumber = numInliers;
			transformations[imageIndex1][imageIndex2].inlierRatio = (double)numInliers / (double)(matches[imageIndex1][imageIndex2-imageIndex1-1].size());
			
			transformations[imageIndex2][imageIndex1].inlierNumber = numInliers;
			transformations[imageIndex2][imageIndex1].inlierRatio = (double)numInliers / (double)(matches[imageIndex1][imageIndex2-imageIndex1-1].size());
		}
		else
		{
			transformations[imageIndex1][imageIndex2].hMatrix = M.getInverseMatrix();
			transformations[imageIndex2][imageIndex1].hMatrix = M;
			
			transformations[imageIndex2][imageIndex1].inlierNumber = numInliers;
			transformations[imageIndex2][imageIndex1].inlierRatio = (double)numInliers / (double)(matches[imageIndex2][imageIndex1-imageIndex2-1].size());
			
			transformations[imageIndex1][imageIndex2].inlierNumber = numInliers;
			transformations[imageIndex1][imageIndex2].inlierRatio = (double)numInliers / (double)(matches[imageIndex2][imageIndex1-imageIndex2-1].size());
		}
	}
	else
	{
		return false;
	}
}


// Compute rigid transforms between all matching images
void GeometryManager::computeTransforms(std::vector<std::vector<Transformation> > &transformations, bool removeBadMatches, vector<vector<keypt_t> > &keyInfo, vector<vector<vector<KeypointMatch> > > &matches) 
{
    unsigned int numImages = matches.size();
	
    for (unsigned int i = 0; i < numImages; i++) 
	{
		for(unsigned int iter=0; iter<matches[i].size(); iter++ )
		{
			bool connect = computeTransform(i, iter+i+1, removeBadMatches, keyInfo, matches, transformations);
			if(!connect)
			{
				if(removeBadMatches)
				{
					matches[i][iter+i+1].clear();
				}
			}
		}
	}
}