/*
 *  KeypointMatcher.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _KEYPOINT_MATCHER_
#define _KEYPOINT_MATCHER_

#include <vector>
#include "Image.h"
#include "Common.h"

using namespace std;

class KeypointMatcher
{
public:
	KeypointMatcher() {}
	std::vector<KeypointMatch> MatchKeys(int num_keys1, vector<short> &k1, 
										 int num_keys2, vector<short> &k2, 
										 double ratio=0.6, int max_pts_visit=200); 
	std::vector<KeypointMatch> CrossMatchKeys(int num_keys1, vector<short> &k1, 
										 int num_keys2, vector<short> &k2, 
										 double ratio=0.6, int max_pts_visit=200); 
	std::vector<KeypointMatch> MatchKeys(int num_keys1, vector<short> &k1, vector<keypt_t> &keyInfo1, 
										 int num_keys2, vector<short> &k2, vector<keypt_t> &keyInfo2, 
										 double ratio, Img &img1, Img &img2) ;
	std::vector<KeypointMatch> MatchKeys(int num_keys1, vector<short> &k1, int num_keys2, vector<short> &k2, 
										 vector<vector<vector<pair<int, int> > > > &featureTracks, vector<vector<int> > &keyStarts, vector<vector<int> > &images,
										 int groupIndex, int maxNoTracks);
	std::vector<KeypointMatch> MatchKeys(int num_keys1, vector<short> &k1, vector<keypt_t> &keyInfo1,
										 int num_keys2, vector<short> &k2, vector<keypt_t> &keyInfo2,
										 double threshold, double &scale, vector<float> &scaleValues, int &noInliers); 
	std::vector<keyPairMatch> MatchKeyPairs(vector<keyPair> &keyPairs1, vector<keyPair> keyPairs2, 
											double ratio=0.6, int max_pts_visit=200);
	std::vector<KeypointMatch> MatchKeysInRegion(int num_keys1, vector<short> &k1, vector<keypt_t> &keyInfo1,
												 int num_keys2, vector<short> &k2, vector<keypt_t> &keyInfo2,
												 vector<Vec2f> centers, int w, int h, double threshold, double &scale, vector<float> &scaleValues) ;
	std::vector<KeypointMatch> MatchKeysOutsideRegion(int num_keys1, vector<short> &k1, vector<keypt_t> &keyInfo1,
													  int num_keys2, vector<short> &k2, vector<keypt_t> &keyInfo2,
													  vector<Vec2f> &centers1, int w1, int h1,
													  vector<Vec2f> &centers2, int w2, int h2) ;
	
	std::vector<KeypointMatch> findUniqueKeys(int num_keys1, vector<short> &k1);
	void findCompatibilityOfMatches(keypt_t pair1a, keypt_t pair1b, keypt_t pair2a, keypt_t pair2b, float &closeness, float &compatibility,
									float closenessThres, float distortionThresh);
	
	void readMatches(const char* filename, vector<vector<vector<KeypointMatch> > > &matches);
	void writeMatches(const char* filename, vector<vector<vector<KeypointMatch> > > &matches);
	
	void setTrackMatches(std::vector<CorrespondenceTrack> &tracks, int numImages, std::vector<std::vector<std::vector<KeypointMatch> > > &trackMatches);
	void computeTracks(vector<vector<keypt_t> > &keyInfo, vector<vector<vector<KeypointMatch> > > &matches, std::vector<CorrespondenceTrack> &tracks);
private:
	//MatchTable m_match;
};


#endif