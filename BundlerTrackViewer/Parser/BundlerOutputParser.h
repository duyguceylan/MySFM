/*
 *  BundlerOutputParser.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include "sfm-driver/sfm.h"
#include "PerspectiveCamera.h"
#include "Image.h"
#include "Common.h"

using namespace std;

class BundlerOutputParser
{
private:
	int noCameras;
	int noTracks;
	
	std::vector<camera_params_t> cameraParams;
	std::vector<CorrespondenceTrack> tracks;
	
public:
	BundlerOutputParser() {}
	~BundlerOutputParser() {} 
	
	void readBundlerOutputFile(const char* filename, int width, int height);
	void writeCamAndImage(const char *outputPath, std::vector<string> &imagePaths, std::vector<PerspectiveCamera> &cameras);
	void readCalibration(const char *outputPath, int noImages, std::vector<PerspectiveCamera> &cameras);
	
	void undistortImage(Img &img, int w, int h, const std::string &out, const camera_params_t &camera);
	
	int getNoOfTrack() {return tracks.size();}
	int getNoOfCorrespondencesInTrack(int index){return tracks[index].correspondences.size(); }
	Correspondence getCorrespondenceInTrack(int trackIndex, int corrIndex) {return tracks[trackIndex].correspondences[corrIndex]; }
	
	void convertSfMOutputToPMVS(const char *input);
};