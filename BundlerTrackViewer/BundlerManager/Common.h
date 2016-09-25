/*
 *  Common.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _COMMON_H
#define _COMMON_H

#include "../MathUtils/MyMatrix.h"

typedef std::pair<int,int> ImageKey;

typedef struct 
{
    double pos[3];
    double color[3];
    std::vector<ImageKey> views;
} point_t;

struct keypt_t{
	bool visited;
    int trackIndex;
	
	int o ;    ///< Keypoint octave index
	
	int ix ;   ///< Keypoint integer X coordinate (unnormalized)
	int iy ;   ///< Keypoint integer Y coordinate (unnormalized)
	int is ;   ///< Keypoint integer scale indiex
	
	float s ;   ///< Keypoint fractional scale index
	
	float x, y;
    float scale;
    float orient;
	bool repeatedKeypt;
	
	bool operator== (const keypt_t &rhs) const
	{
		if(fabs(x-rhs.x)<0.5 && fabs(y-rhs.y)<0.5 && fabs(orient-rhs.orient)<0.001)
			return true;
		else
			return false;
	}
	
};

struct keyPair
{
	keypt_t keyInfo1;
	keypt_t keyInfo2;
	int keyIdx1;
	int keyIdx2;
	float angle;
	std::vector<short> descriptor;
	int descriptorComputed;
	keyPair(keypt_t key1, keypt_t key2, int idx1, int idx2) : keyInfo1(key1), keyInfo2(key2), keyIdx1(idx1), keyIdx2(idx2) {descriptorComputed = -1; descriptor.resize(256);}
	keyPair() {}
	
};

struct keyPairMatch
{
	int idx1;
	int idx2;
	float confidence;
	
	keyPairMatch(int keyPair1, int keyPair2) : idx1(keyPair1), idx2(keyPair2) {}
	keyPairMatch() {}
	
	bool operator<(const keyPairMatch &rhs) const
	{
		if(confidence < rhs.confidence)
			return true;
		else
			return false;
	}
};

struct Correspondence
{
	int imageIndex;
	int keyIndex;
	Vec2f position;
};

struct CorrespondenceTrack
{
	Vec3f pos3D;
	Vec3uc color;
	vector<Correspondence> correspondences;
};

typedef struct Transformation 
{
	Matrix3f fMatrix;
	
	Matrix3f hMatrix;
	int inlierNumber;
	double inlierRatio;
	
	/* double m_ematrix[9];
	 //color correction
	 double m_gain[3];
	 double m_bias[3];
	 */
	bool valid;	
};

/* Data struct for matches */
class KeypointMatch {
public:
    KeypointMatch() { validMatch = true; };
	
	KeypointMatch(int idx1, int idx2) : m_idx1(idx1), m_idx2(idx2) { validMatch = true;}
	
	bool operator==(const KeypointMatch &rhs) const
	{
		if(rhs.m_idx1 == m_idx1 && rhs.m_idx2 == m_idx2)
			return true;
		else
			return false;
	}
	
	int m_idx1, m_idx2;
	bool validMatch;
	bool repeatedMatch;
};


#endif