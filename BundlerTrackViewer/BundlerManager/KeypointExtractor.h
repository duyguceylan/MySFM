/*
 *  KeyPointExtractor.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/16/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include <zlib.h>
#include "Common.h"
#include "SIFT.h"

using namespace std;

class KeypointExtractor
{
public:
	KeypointExtractor() {}
	void WriteKeyFile(const char *filename, vector<short> &keys, vector<keypt_t> &info);
	int ReadKeyFile(const char *filename, vector<short> &keys, vector<keypt_t> &info);
	int ReadKeysGzip(gzFile fp, vector<short> &keys, vector<keypt_t> &info);
	int ReadKeys(FILE *fp, vector<short> &keys, vector<keypt_t> &info);
	
	void formFeaturePairs(std::vector<keypt_t> &keyInfo, std::vector<short> &keys, std::vector<keyPair> &keyPairs, float minDist, float maxDist);
};