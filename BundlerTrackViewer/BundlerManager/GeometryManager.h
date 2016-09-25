/*
 *  GeometryManager.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _GEOMETRY_MANAGER_H
#define _GEOMETRY_MANAGER_H

#include <vector>
#include "EpipolarGeometry.h"
#include "Common.h"

using namespace std;

class GeometryManager
{
public:
	GeometryManager();
	
	bool computeEpipolarGeometry(int imageIndex1, int imageIndex2, bool removeBadMatches,vector<vector<keypt_t> > &keyInfo, 
								 vector<vector<vector<KeypointMatch> > > &matches, std::vector<std::vector<Transformation> > &transformations);
	void computeEpipolarGeometry(std::vector<std::vector<Transformation> > &transformations, bool removeBadMatches, vector<vector<keypt_t> > &keyInfo, 
								 vector<vector<vector<KeypointMatch> > > &matches);
	
	bool computeTransform(int imageIndex1, int imageIndex2, bool removeBadMatches,vector<vector<keypt_t> > &keyInfo, 
								 vector<vector<vector<KeypointMatch> > > &matches, std::vector<std::vector<Transformation> > &transformations);
	void computeTransforms(std::vector<std::vector<Transformation> > &transformations, bool removeBadMatches, vector<vector<keypt_t> > &keyInfo, 
								 vector<vector<vector<KeypointMatch> > > &matches);

private:
	EpipolarGeometry *eg;
	
	int fMatrixMaxTrialNumber;
	int fMatrixMinInlierNumber;
	double fMatrixThreshold;
	
	int homographyTrialNumber;
	int homographyInlierNumber;
	double homographyThreshold;
};

#endif