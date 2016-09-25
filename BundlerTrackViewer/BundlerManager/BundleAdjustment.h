/*
 *  BundleAdjustment.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _BUNDLE_ADJUSTMENT_H
#define _BUNDLE_ADJUSTMENT_H

#include <vector>
#include "sfm-driver/sfm.h"
#include "Common.h"

class BundleAdjustment
{
public:
	BundleAdjustment();
	void setCameraConstraints(int imgIndex, camera_params_t &cam);
	void initializeCameraParams(camera_params_t &camera);
	void setupInitialCameraPair(int imageIndex1, int imageIndex2, double &initFocalLength1, double &initFocalLength2,
								std::vector<camera_params_t> &cameras, std::vector<std::vector<keypt_t> > &keyInfo, 
								std::vector<std::vector<std::vector<KeypointMatch> > > &trackMatches, 
								std::vector<CorrespondenceTrack> &tracks);
	
	void pickInitialPair(int &imgIndex1, int &imgIndex2, bool useInitFocalOnly, std::vector<std::vector<std::vector<KeypointMatch> > > &trackMatches,
						 std::vector<std::vector<Transformation> > &transformations);
private:
	double projectionEstimationThreshold;
};

#endif