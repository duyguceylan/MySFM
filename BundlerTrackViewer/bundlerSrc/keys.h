/* 
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* keys.h */
/* Class for SIFT keypoints */

#ifndef __keys_h__
#define __keys_h__

#include <vector>
#include <stdio.h>
#include <zlib.h>
#include "Common.h"

#include "ANN.h"

// #define KEY_LIMIT 65536

class Keypoint {
public:    
    Keypoint()  
    { m_x = 0.0; m_y = 0.0; m_extra = -1; m_track = -1; }

    Keypoint(float x, float y) :
	m_x(x), m_y(y)
    { m_r = 0; m_g = 0; m_b = 0; }

    virtual ~Keypoint() {}   
    
    virtual unsigned char *GetDesc() {
        return NULL;
    }

    float m_x, m_y;        /* Subpixel location of keypoint. */
    unsigned char m_r, m_g, m_b;          /* Color of this key */

    int m_extra;  /* 4 bytes of extra storage */
    int m_track;  /* Track index this point corresponds to */
	bool repeatedKeypt;
};

class KeypointWithDesc : public Keypoint
{
public:
    virtual ~KeypointWithDesc() {}
    
    virtual unsigned char *GetDesc() {
        return m_d;
    }

    KeypointWithDesc()  
    { m_x = 0.0; m_y = 0.0; m_d = NULL; 
        m_extra = -1; m_track = -1; }

    KeypointWithDesc(float x, float y, unsigned char *d) :
	Keypoint(x, y), m_d(d)
    { }

    unsigned char *m_d;    /* Vector of descriptor values */
};

class KeypointWithScaleRot : public KeypointWithDesc
{
public:
    virtual ~KeypointWithScaleRot() {}
	
	KeypointWithScaleRot()  
    { m_x = 0.0; m_y = 0.0; m_d = NULL; 
        m_extra = -1; m_track = -1; m_scale = 0.0; m_orient = 0.0;}
	
	KeypointWithScaleRot(float x, float y, float scale, float ori, unsigned char *d) :
	KeypointWithDesc(x,y,d), m_scale(scale), m_orient(ori)
    { }

    float m_scale, m_orient;
};

/* Data struct for matches (with score field) */
class KeypointMatchWithScore : public KeypointMatch {
public:
    KeypointMatchWithScore() 
    { }

#ifdef KEY_LIMIT
    KeypointMatchWithScore (unsigned short idx1, unsigned short idx2, 
                            float score) :
#else
    KeypointMatchWithScore (int idx1, int idx2, float score) :
#endif
	KeypointMatch(idx1, idx2), m_score(score)
    { }

#ifdef KEY_LIMIT
    unsigned short m_idx1, m_idx2;
#else
    int m_idx1, m_idx2;
#endif

    float m_score;
};

/* Returns the number of keys in a file */
int GetNumberOfKeys(const char *filename);

/* This reads a keypoint file from a given filename and returns the list
 * of keypoints. */
std::vector<Keypoint> ReadKeyFile(const char *filename);
std::vector<KeypointWithDesc> 
    ReadKeyFileWithDesc(const char *filename, bool descriptor);
std::vector<KeypointWithScaleRot> 
    ReadKeyFileWithScaleRot(const char *filename, bool descriptor);
int ReadKeyFile(const char *filename, short **keys, 
                keypt_t **info = NULL);
int ReadKeyPositions(const char *filename, keypt_t **info);

/* Read keypoints from the given file pointer and return the list of
 * keypoints.  The file format starts with 2 integers giving the total
 * number of keypoints and the size of descriptor vector for each
 * keypoint (currently assumed to be 128). Then each keypoint is
 * specified by 4 floating point numbers giving subpixel row and
 * column location, scale, and orientation (in radians from -PI to
 * PI).  Then the descriptor vector for each keypoint is given as a
 * list of integers in range [0,255]. */
std::vector<Keypoint> ReadKeys(FILE *fp, bool descriptor);
int ReadKeys(FILE *fp, short **keys, keypt_t **info = NULL);

/* Read keys more quickly */
std::vector<KeypointWithDesc> ReadKeysFast(FILE *fp, bool descriptor, 
                                           float **scale = NULL, 
                                           float **orient = NULL,
										   bool **repeated = NULL);
std::vector<KeypointWithDesc> ReadKeysFastGzip(gzFile fp, bool descriptor,
                                               float **scale = NULL, 
                                               float **orient = NULL);
int ReadKeysGzip(gzFile fp, short **keys, keypt_t **info = NULL);

/* Read keys from binary file */
std::vector<KeypointWithDesc> ReadKeysFastBin(FILE *fp, bool descriptor,
                                                 float **scales = NULL, 
                                                 float **orients = NULL);

/* Read keys from gzipped binary file */
std::vector<KeypointWithDesc> ReadKeysFastBinGzip(gzFile fp, bool descriptor,
                                                  float **scales = NULL, 
                                                  float **orients = NULL);
/* Read keys using MMAP to speed things up */
std::vector<Keypoint> ReadKeysMMAP(FILE *fp);

/* Create a search tree for the given set of keypoints */
ANNkd_tree *CreateSearchTree(int num_keys, unsigned char *keys);

/* Compute likely matches between two sets of keypoints */
std::vector<KeypointMatch> MatchKeys(const std::vector<KeypointWithDesc> &k1, 
				     const std::vector<KeypointWithDesc> &k2,
				     bool registered = false, 
				     double ratio = 0.6);
std::vector<KeypointMatch> MatchKeys(int num_keys1, short *k1, 
									 int num_keys2, short *k2,
									 double ratio = 0.6, 
                                     int max_pts_visit = 200);

std::vector<KeypointMatch> MatchKeys(int num_keys1, unsigned char *k1, 
                                     ANNkd_tree *tree2,
									 double ratio = 0.6, 
                                     int max_pts_visit = 200);

/* Compute likely matches between two sets of keypoints */
std::vector<KeypointMatchWithScore> 
    MatchKeysWithScore(const std::vector<KeypointWithDesc> &k1, 
                       const std::vector<KeypointWithDesc> &k2,
                       bool registered = false, 
                       double ratio = 0.6);

/* Prune matches so that they are 1:1 */
std::vector<KeypointMatchWithScore> 
    PruneMatchesWithScore(const std::vector<KeypointMatchWithScore> &matches);

std::vector<KeypointMatch> 
    MatchKeysExhaustive(const std::vector<KeypointWithDesc> &k1, 
                        const std::vector<KeypointWithDesc> &k2,
                        bool registered = false, 
                        double ratio = 0.6);

#endif /* __keys_h__ */
