/*
 *  KeypointMatcher.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include <queue>
#include "KeypointMatcher.h"
#include "ANN.h"

bool compareFirst(const KeypointMatch &k1, const KeypointMatch &k2) 
{
    return (k1.m_idx1 < k2.m_idx1);
}

std::vector<KeypointMatch> KeypointMatcher::MatchKeys(int num_keys1, vector<short> &k1, vector<keypt_t> &keyInfo1, int num_keys2, vector<short> &k2, vector<keypt_t> &keyInfo2, 
													  double ratio, Img &img1, Img &img2) 
{
	std::vector<KeypointMatch> matches;
	
	ANNpointArray pts = annAllocPts(num_keys2, 128);
	for (int k = 0; k < num_keys2; k++) 
	{
		Color c = img2(keyInfo2[k].x, keyInfo2[k].y);
		if((c-Color(1.0,1.0,1.0)).length() > 0.02)
			continue;
		for(int j=0; j<128; j++)
				pts[k][j] = k2[128*k+j];
	}
	
	ANNkd_tree *tree = new ANNkd_tree(pts, num_keys2, 128, 16);
	
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
		Color c = img1(keyInfo1[i].x, keyInfo1[i].y);
		
		if((c-Color(1.0,1.0,1.0)).length() > 0.02)
			continue;
		
		//if((c-Color(0.0,0.0,0.0)).length() < 0.02)
		//	continue;
		
		// Create a new array of points
		/*ANNpointArray pts = annAllocPts(num_keys2, 128);
		vector<int> indices;
		int count = 0;
		for (int k = 0; k < num_keys2; k++) 
		{
			Color nc = img2(keyInfo2[k].x, keyInfo2[k].y);
			
			if((c-nc).length() < 0.02)
			{
				for(int j=0; j<128; j++)
					pts[count][j] = k2[128*k+j];
				indices.push_back(k);
				count++;
			}
			
		}
		
		if(count < 2)
		{
			annDeallocPts(pts);
			continue;
		}
		
		// Create a search tree for k2
		ANNkd_tree *tree = new ANNkd_tree(pts, count, 128, 16);
		*/
		
        ANNidxArray nn_idx = new ANNidx[2];
        ANNdistArray dist = new ANNdist[2];
		ANNpoint pt = annAllocPt(128);
		for(int j=0; j<128; j++)
			pt[j] = k1[128*i+j];
		
		tree->annkSearch(pt, 2, nn_idx, dist, 0.0);
		
		if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
		{
			Color nc = img2(keyInfo2[nn_idx[0]].x, keyInfo2[nn_idx[0]].y);
			
			if((c-nc).length() < 0.02)
			{
				
				matches.push_back(KeypointMatch(i, /*indices[*/nn_idx[0]/*]*/));
			}
        }
		
		delete [] nn_idx;
		delete [] dist;
		annDeallocPt(pt);
		/* Cleanup */
		//annDeallocPts(pts);
		
		//delete tree;
    }
	
	printf("%d matches found\n", matches.size());
	
	annDeallocPts(pts);
	
	delete tree;
	
    return matches;
}

std::vector<KeypointMatch> KeypointMatcher::MatchKeys(int num_keys1, vector<short> &k1, int num_keys2, vector<short> &k2, 
													  double ratio, int max_pts_visit) 
{
	annMaxPtsVisit(max_pts_visit);
	
    std::vector<KeypointMatch> matches;
	
	if(num_keys1 == 0 || num_keys2 == 0)
		return matches;
	
    // Create a new array of points
    ANNpointArray pts = annAllocPts(num_keys2, 128);
	
    for (int i = 0; i < num_keys2; i++) 
	{
		for(int j=0; j<128; j++)
			pts[i][j] = k2[128*i+j];
    }
	
    // Create a search tree for k2
    ANNkd_tree *tree = new ANNkd_tree(pts, num_keys2, 128, 16);
    
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
        ANNidxArray nn_idx = new ANNidx[2];
        ANNdistArray dist = new ANNdist[2];
		ANNpoint pt = annAllocPt(128);
		for(int j=0; j<128; j++)
			pt[j] = k1[128*i+j];
		
		tree->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
		
		if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
		{
			matches.push_back(KeypointMatch(i, nn_idx[0]));
        }
		
		delete [] nn_idx;
		delete [] dist;
		annDeallocPt(pt);
    }
	
	printf("%d matches found\n", matches.size());
	
    /* Cleanup */
    annDeallocPts(pts);
	
    delete tree;
	
    return matches;
}

std::vector<KeypointMatch> KeypointMatcher::MatchKeys(int num_keys1, vector<short> &k1, int num_keys2, vector<short> &k2, 
													  vector<vector<vector<pair<int, int> > > > &featureTracks, vector<vector<int> > &keyStarts, vector<vector<int> > &images,
													  int groupIndex, int maxNoTracks)
{
	double ratio = 0.6;
	double distThreshold = 3000;
	int max_pts_visit = 200;
	
	annMaxPtsVisit(max_pts_visit);
	
    std::vector<KeypointMatch> matches;
	
    // Create a new array of points
    ANNpointArray pts = annAllocPts(num_keys2, 128);
	
    for (int i = 0; i < num_keys2; i++) 
	{
		for(int j=0; j<128; j++)
			pts[i][j] = k2[128*i+j];
    }
	
    // Create a search tree for k2
    ANNkd_tree *tree = new ANNkd_tree(pts, num_keys2, 128, 16);
    
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
		ANNpoint pt = annAllocPt(128);
		for(int j=0; j<128; j++)
			pt[j] = k1[128*i+j];
		
		if(groupIndex == -1)
		{
			ANNidxArray nn_idx = new ANNidx[2];
			ANNdistArray dist = new ANNdist[2];
			tree->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
			
			if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
			{
				matches.push_back(KeypointMatch(i, nn_idx[0]));
			}
			delete [] nn_idx;
			delete [] dist;
		}
		else
		{
			maxNoTracks = 2;
			
			ANNidxArray nn_idx = new ANNidx[maxNoTracks];
			ANNdistArray dist = new ANNdist[maxNoTracks];
			tree->annkPriSearch(pt, maxNoTracks, nn_idx, dist, 0.0);
			
			if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
			{
				matches.push_back(KeypointMatch(i, nn_idx[0]));
			}
			else
			{
				int img1;
				int feature1;
				for(int k=0; k<keyStarts[groupIndex].size(); k++)
				{
					if(nn_idx[0] < keyStarts[groupIndex][k])
					{
						img1 = images[groupIndex][k];
						feature1 = nn_idx[0];
						if(k<0)
							feature1 = nn_idx[0] - keyStarts[groupIndex][k-1]; 
						break;
					}
				}
				
				int img2;
				int feature2;
				for(int k=0; k<keyStarts[groupIndex].size(); k++)
				{
					if(nn_idx[1] < keyStarts[groupIndex][k])
					{
						img2 = images[groupIndex][k];
						feature2 = nn_idx[1];
						if(k<0)
							feature2 = nn_idx[1] - keyStarts[groupIndex][k-1]; 
						break;
					}
				}
				
				bool found = false;
				for(int t=0; t<featureTracks[groupIndex].size(); t++)
				{
					for(int k=0; k<featureTracks[groupIndex][t].size(); k++)
					{
						if(featureTracks[groupIndex][t][k].first == img1 && featureTracks[groupIndex][t][k].second == feature1)
						{
							for(int e1=0; e1<featureTracks[groupIndex][t].size(); e1++)
							{
								if(featureTracks[groupIndex][t][k].first == img2 && featureTracks[groupIndex][t][k].second == feature2)
								{
									found = true;
									break;
								}
							}
						}
						if(found)
							break;
					}
					if(found)
						break;
				}
				if (found) 
				{
					matches.push_back(KeypointMatch(i, nn_idx[0]));
				}
				delete [] nn_idx;
				delete [] dist;
				
				annDeallocPt(pt);
			}
		}
	}
			
			/*int matchIndex = -1;
			float minDist = 0;
			bool trackFound = false;
			bool add=true;
			float avgDist = 0.0;
			for(int m=0; m<1; m++)
			{
				//find image index
				int img;
				int feature;
				for(int k=0; k<keyStarts[groupIndex].size(); k++)
				{
					if(nn_idx[0] < keyStarts[groupIndex][k])
					{
						img = images[groupIndex][k];
						feature = nn_idx[m];
						if(k<0)
							feature = nn_idx[m] - keyStarts[groupIndex][k-1]; 
						break;
					}
				}
				
				
				for(int t=0; t<featureTracks[groupIndex].size(); t++)
				{
					for(int k=0; k<featureTracks[groupIndex][t].size(); k++)
					{
						if(featureTracks[groupIndex][t][k].first == img && featureTracks[groupIndex][t][k].second == feature)
						{
							
							for(int e1=0; e1<featureTracks[groupIndex][t].size(); e1++)
							{
								int img2 = featureTracks[groupIndex][t][e1].first;
								int f2 = featureTracks[groupIndex][t][e1].second;
								img2 = find(images[groupIndex].begin(), images[groupIndex].end(), img2) - images[groupIndex].begin();
								if(img2 > 0)
									f2 = f2 + keyStarts[groupIndex][img2-1];
								float tmpDist = 0.0;
								for(int d=0; d<128; d++)
								{
									tmpDist += (k1[128*i+d] - k2[128*f2+d])*(k1[128*i+d] - k2[128*f2+d]);
								}
								if(tmpDist > distThreshold)
								{
									add= false;
									break;
								}
								avgDist += tmpDist;
							}
							avgDist /= (float)(featureTracks[groupIndex][t].size());
							trackFound = true;
							break;
						}
					}
					if(trackFound)
						break;
				}
			}
			if(trackFound)		
			{
				if(add)
				{
					matches.push_back(KeypointMatch(i, nn_idx[0]));
				}
			}
			else if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
			{
				matches.push_back(KeypointMatch(i, nn_idx[0]));
			}
			
			delete [] nn_idx;
			delete [] dist;
		}
				
		annDeallocPt(pt);
    }
	
	printf("%d matches found\n", matches.size());
	
    /* Cleanup */
    annDeallocPts(pts);
	
    delete tree;
	
    return matches;
}


std::vector<KeypointMatch> KeypointMatcher::MatchKeys(int num_keys1, vector<short> &k1, vector<keypt_t> &keyInfo1,
													  int num_keys2, vector<short> &k2, vector<keypt_t> &keyInfo2,
													  double threshold, double &scale, vector<float> &scaleValues, int &noInliers) 
{
	int max_pts_visit = 200;
	float ratio = 0.6;
	annMaxPtsVisit(max_pts_visit);
	
    std::vector<KeypointMatch> matches;
	
	if(num_keys1 == 0 || num_keys2 == 0)
		return matches;
	
    // Create a new array of points
    ANNpointArray pts = annAllocPts(num_keys2, 128);
	
    for (int i = 0; i < num_keys2; i++) 
	{
		for(int j=0; j<128; j++)
			pts[i][j] = k2[128*i+j];
    }
	
	float totalScale = 0.0;
	
    // Create a search tree for k2
    ANNkd_tree *tree = new ANNkd_tree(pts, num_keys2, 128, 16);
    
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
		ANNidxArray nn_idx = new ANNidx[2];
		ANNdistArray dist = new ANNdist[2];
		ANNpoint pt = annAllocPt(128);
		for(int j=0; j<128; j++)
			pt[j] = k1[128*i+j];
				
		tree->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
				
		//printf("ratio:%f\n", (double) dist[0] / (double) dist[1]);
		if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
		{
			totalScale += keyInfo2[nn_idx[0]].scale/keyInfo1[i].scale;
			scaleValues.push_back(keyInfo2[nn_idx[0]].scale/keyInfo1[i].scale);
			matches.push_back(KeypointMatch(i, nn_idx[0]));
		}
				
		delete [] nn_idx;
		delete [] dist;
		annDeallocPt(pt);
	}
	
	sort(scaleValues.begin(), scaleValues.end());
	vector<vector<float> > clusters;
	float lastValue;
	int maxClusterSize = -1;
	int maxClusterIndex = -1;
	for(int loop=0; loop<scaleValues.size(); loop++)
	{
		lastValue = scaleValues[loop];
		//printf("scale:%f\n", lastValue);
		int count = 0;
		for(int i=0; i<scaleValues.size(); i++)
		{
			//printf("s:%f\n", scaleValues[i]);
			if(abs(scaleValues[i] - lastValue) <= threshold)
			{
				count++;
			}
		}
		if(maxClusterSize == -1)
		{
			maxClusterSize = count;
			maxClusterIndex = loop;
		}
		else if(count > maxClusterSize)
		{
			maxClusterSize = count;
			maxClusterIndex = loop;
		}
	}
	int thresholdClusterSize = 40;
	noInliers = 0;
	float scaleRatio = (float)maxClusterSize / (float)scaleValues.size();
	if(maxClusterIndex != -1 /*&& scaleRatio > 0.5maxClusterSize >= thresholdClusterSize*/)
	{
		lastValue = scaleValues[maxClusterIndex];
		vector<float> cluster;
		double avgScale = 0.0;
		for(int i=0; i<scaleValues.size(); i++)
		{
			if(abs(scaleValues[i] - lastValue) <= threshold)
			{
				cluster.push_back(scaleValues[i]);
				avgScale += scaleValues[i];
			}
		}
		sort(cluster.begin(), cluster.end());
		scale = cluster[cluster.size()/2];//clusters[maxClusterIndex][maxClusterSize/2];//scaleValues[scaleValues.size()/2];
		printf("scale:%f\n", scale);
		noInliers = maxClusterSize;
	}
	else
		scale = 0.0;
	//printf("%d matches found\n", matches.size());
	
    // Cleanup 
    annDeallocPts(pts);
	
    delete tree;
	
	return matches;
}


std::vector<KeypointMatch> KeypointMatcher::CrossMatchKeys(int num_keys1, vector<short> &k1, int num_keys2, vector<short> &k2, 
														   double ratio, int max_pts_visit) 
{
	annMaxPtsVisit(max_pts_visit);
	
    std::vector<KeypointMatch> matches;
	
	// Create a new array of points
    ANNpointArray pts1 = annAllocPts(num_keys1, 128);
	
    for (int i = 0; i < num_keys1; i++) 
	{
		for(int j=0; j<128; j++)
			pts1[i][j] = k1[128*i+j];
    }
	
    // Create a search tree for k2
    ANNkd_tree *tree1 = new ANNkd_tree(pts1, num_keys1, 128, 16);
	
    // Create a new array of points
    ANNpointArray pts2 = annAllocPts(num_keys2, 128);
	
    for (int i = 0; i < num_keys2; i++) 
	{
		for(int j=0; j<128; j++)
			pts2[i][j] = k2[128*i+j];
    }
	
    // Create a search tree for k2
    ANNkd_tree *tree2 = new ANNkd_tree(pts2, num_keys2, 128, 16);
    
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
        ANNidxArray nn_idx = new ANNidx[2];
        ANNdistArray dist = new ANNdist[2];
		ANNpoint pt = annAllocPt(128);
		for(int j=0; j<128; j++)
			pt[j] = k1[128*i+j];
		
		tree2->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
		
		if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
		{
			ANNidxArray nn_idx2 = new ANNidx[1];
			ANNdistArray dist2 = new ANNdist[1];
			ANNpoint pt2 = annAllocPt(128);
			for(int j=0; j<128; j++)
				pt2[j] = k2[128*nn_idx[0]+j];
			
			tree1->annkPriSearch(pt2, 1, nn_idx2, dist2, 0.0);
			if(nn_idx2[0] == i)
				matches.push_back(KeypointMatch(i, nn_idx[0]));
			
			delete [] nn_idx2;
			delete [] dist2;
			annDeallocPt(pt2);
        }
		
		delete [] nn_idx;
		delete [] dist;
		annDeallocPt(pt);
    }
	
	printf("%d matches found\n", matches.size());
	
    /* Cleanup */
    annDeallocPts(pts1);
    delete tree1;
	
	annDeallocPts(pts2);
    delete tree2;
	
    return matches;
}

std::vector<keyPairMatch> KeypointMatcher::MatchKeyPairs(vector<keyPair> &keyPairs1, vector<keyPair> keyPairs2, 
														  double ratio, int max_pts_visit) 
{
	annMaxPtsVisit(max_pts_visit);
	
    std::vector<keyPairMatch> matches;
	
    // Create a new array of points
    ANNpointArray pts = annAllocPts(keyPairs2.size(), 256);
	
    for (int i = 0; i < keyPairs2.size(); i++) 
	{
		for(int j=0; j<256; j++)
			pts[i][j] = keyPairs2[i].descriptor[j];
    }
	
    // Create a search tree for k2
    ANNkd_tree *tree = new ANNkd_tree(pts, keyPairs2.size(), 256, 16);
    
	// Now do the search
    for (int i = 0; i < keyPairs1.size(); i++) 
	{
        ANNidxArray nn_idx = new ANNidx[2];
        ANNdistArray dist = new ANNdist[2];
		ANNpoint pt = annAllocPt(256);
		for(int j=0; j<256; j++)
		{
			pt[j] = keyPairs1[i].descriptor[j];
			//printf("%f\n", pt[j]);
		}
		
		tree->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
		
		if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
		{
			keyPairMatch m(i, nn_idx[0]);
			m.confidence = ((double) dist[0]) / ((double) dist[1]);
			matches.push_back(m);
        }
		
		delete [] nn_idx;
		delete [] dist;
		annDeallocPt(pt);
    }
	
	printf("%d matches found\n", matches.size());
	
    /* Cleanup */
    annDeallocPts(pts);
	
    delete tree;
	
    return matches;
}


std::vector<KeypointMatch> KeypointMatcher::MatchKeysInRegion(int num_keys1, vector<short> &k1, vector<keypt_t> &keyInfo1,
															  int num_keys2, vector<short> &k2, vector<keypt_t> &keyInfo2,
															  vector<Vec2f> centers, int w, int h, double threshold, double &scale,
															  vector<float> &scaleValues) 
{
	int max_pts_visit = 200;
	float ratio = 0.6;
	annMaxPtsVisit(max_pts_visit);
	
    std::vector<KeypointMatch> matches;
	
    // Create a new array of points
    ANNpointArray pts = annAllocPts(num_keys2, 128);
	
    for (int i = 0; i < num_keys2; i++) 
	{
		for(int j=0; j<128; j++)
			pts[i][j] = k2[128*i+j];
    }
	
	float totalScale = 0.0;
	
    // Create a search tree for k2
    ANNkd_tree *tree = new ANNkd_tree(pts, num_keys2, 128, 16);
    
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
		float x = keyInfo1[i].x;
		float y = keyInfo1[i].y;
		for(int c=0; c<centers.size(); c++)
		{
			if(x>=centers[c][0]-w/2 && x<centers[c][0]+w/2 && y>=centers[c][1]-h/2 && y<centers[c][1]+h/2)
			{
				ANNidxArray nn_idx = new ANNidx[2];
				ANNdistArray dist = new ANNdist[2];
				ANNpoint pt = annAllocPt(128);
				for(int j=0; j<128; j++)
					pt[j] = k1[128*i+j];
		
				tree->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
		
				if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
				{
					totalScale += keyInfo2[nn_idx[0]].scale/keyInfo1[i].scale;
					scaleValues.push_back(keyInfo2[nn_idx[0]].scale/keyInfo1[i].scale);
					matches.push_back(KeypointMatch(i, nn_idx[0]));
				}
		
				delete [] nn_idx;
				delete [] dist;
				annDeallocPt(pt);
				break;
			}
		}
    }
	
	sort(scaleValues.begin(), scaleValues.end());
	vector<vector<float> > clusters;
	float lastValue;
	int maxClusterSize = -1;
	int maxClusterIndex = -1;
	for(int loop=0; loop<scaleValues.size(); loop++)
	{
		lastValue = scaleValues[loop];
		printf("scale:%f\n", lastValue);
		int count = 0;
		for(int i=0; i<scaleValues.size(); i++)
		{
			//printf("s:%f\n", scaleValues[i]);
			if(abs(scaleValues[i] - lastValue) <= threshold)
			{
				count++;
			}
		}
		if(maxClusterSize == -1)
		{
			maxClusterSize = count;
			maxClusterIndex = loop;
		}
		else if(count > maxClusterSize)
		{
			maxClusterSize = count;
			maxClusterIndex = loop;
		}
	}
	if(maxClusterIndex != -1 && maxClusterSize >= 10)
	{
		lastValue = scaleValues[maxClusterIndex];
		vector<float> cluster;
		for(int i=0; i<scaleValues.size(); i++)
		{
			if(scaleValues[i] - lastValue <= threshold)
			{
				cluster.push_back(scaleValues[i]);
			}
		}
		sort(cluster.begin(), cluster.end());
		scale = cluster[cluster.size()/2];//clusters[maxClusterIndex][maxClusterSize/2];//scaleValues[scaleValues.size()/2];
		printf("scale:%f\n", scale);
	}
	else
		scale = 0.0;
	//printf("%d matches found\n", matches.size());
	
    // Cleanup 
    annDeallocPts(pts);
	
    delete tree;

	return matches;
}

std::vector<KeypointMatch> KeypointMatcher::MatchKeysOutsideRegion(int num_keys1, vector<short> &k1, vector<keypt_t> &keyInfo1,
																   int num_keys2, vector<short> &k2, vector<keypt_t> &keyInfo2,
																   vector<Vec2f> &centers1, int w1, int h1,
																   vector<Vec2f> &centers2, int w2, int h2) 
{
	int max_pts_visit = 200;
	float ratio = 0.6;
	annMaxPtsVisit(max_pts_visit);
	
    std::vector<KeypointMatch> matches;
	
	
    // Create a new array of points
    ANNpointArray pts = annAllocPts(num_keys2, 128);
	
    for (int i = 0; i < num_keys2; i++) 
	{
		for(int j=0; j<128; j++)
			pts[i][j] = k2[128*i+j];
    }
	
    // Create a search tree for k2
    ANNkd_tree *tree = new ANNkd_tree(pts, num_keys2, 128, 16);
    
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
		float x = keyInfo1[i].x;
		float y = keyInfo1[i].y;
		bool valid = true;
		for(int c=0; c<centers1.size(); c++)
		{
			if(x>=centers1[c][0]-w1/2 && x<centers1[c][0]+w1/2 && y>=centers1[c][1]-h1/2 && y<centers1[c][1]+h1/2)
			{
				valid = false;
				//printf("not valid\n");
				break;
			}
		}
		if(valid)	
		{
			ANNidxArray nn_idx = new ANNidx[2];
			ANNdistArray dist = new ANNdist[2];
			ANNpoint pt = annAllocPt(128);
			for(int j=0; j<128; j++)
				pt[j] = k1[128*i+j];
			
			tree->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
			
			if (((double) dist[0]) < ratio * ratio * ((double) dist[1])) 
			{
				x = keyInfo2[nn_idx[0]].x;
				y = keyInfo2[nn_idx[0]].y;
				for(int c=0; c<centers2.size(); c++)
				{
					if(x>=centers2[c][0]-w2/2 && x<centers2[c][0]+w2/2 && y>=centers2[c][1]-h2/2 && y<centers2[c][1]+h2/2)
					{
						valid = false;
						//printf("not valid\n");
						break;
					}
				}
				
				if(valid)
				{
					matches.push_back(KeypointMatch(i, nn_idx[0]));
				}
			}
				
			delete [] nn_idx;
			delete [] dist;
			annDeallocPt(pt);
		}
    }
	
	 // Cleanup 
    annDeallocPts(pts);
	
    delete tree;
	
	return matches;
}

std::vector<KeypointMatch> KeypointMatcher::findUniqueKeys(int num_keys1, vector<short> &k1) 
{
	int max_pts_visit = 200;
	
	annMaxPtsVisit(max_pts_visit);
	
    std::vector<KeypointMatch> matches;
	
	double avgDistance = 0.0;
	
    // Create a new array of points
    ANNpointArray pts = annAllocPts(num_keys1, 128);
	
    for (int i = 0; i < num_keys1; i++) 
	{
		for(int j=0; j<128; j++)
			pts[i][j] = k1[128*i+j];
    }
	
    // Create a search tree for k2
    ANNkd_tree *tree = new ANNkd_tree(pts, num_keys1, 128, 16);
    
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
		ANNidxArray nn_idx = new ANNidx[2];
		ANNdistArray dist = new ANNdist[2];
		ANNpoint pt = annAllocPt(128);
		for(int j=0; j<128; j++)
			pt[j] = k1[128*i+j];
			
		tree->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
		
		avgDistance += (double) dist[1];
	
		delete [] nn_idx;
		delete [] dist;
		annDeallocPt(pt);
    }
	
	avgDistance /= num_keys1;
	double thresholdDistance = avgDistance;
	
	// Now do the search
    for (int i = 0; i < num_keys1; i++) 
	{
		ANNidxArray nn_idx = new ANNidx[2];
		ANNdistArray dist = new ANNdist[2];
		ANNpoint pt = annAllocPt(128);
		for(int j=0; j<128; j++)
			pt[j] = k1[128*i+j];
		
		tree->annkPriSearch(pt, 2, nn_idx, dist, 0.0);
		
		if((double)(dist[1]) > thresholdDistance)
			matches.push_back(KeypointMatch(i, i));
		
		delete [] nn_idx;
		delete [] dist;
		annDeallocPt(pt);
    }
	
	// Cleanup 
    annDeallocPts(pts);
	
    delete tree;
	
	return matches;
}

void KeypointMatcher::findCompatibilityOfMatches(keypt_t pair1a, keypt_t pair1b, keypt_t pair2a, keypt_t pair2b, float &closeness, float &compatibility,
												 float closenessThres, float distortionThresh)
{
	float dist1 = (Vec2f(pair1a.x, pair1a.y)-Vec2f(pair2a.x, pair2a.y)).length();
	float dist2 = (Vec2f(pair1b.x, pair1b.y)-Vec2f(pair2b.x, pair2b.y)).length();
	if(dist1<closenessThres || dist2<closenessThres)
		closeness = 1.0;
	else
		closeness = 0.0;
	
	if(fabs(dist1-dist2) < distortionThresh)
		compatibility = 1.0;
	else
		compatibility = 0.0;
}

void KeypointMatcher::writeMatches(const char* filename, vector<vector<vector<KeypointMatch> > > &matches)
{
	FILE *f = fopen(filename, "w");
	
	int numImages= (int) matches.size();
	
	for(int i=0; i<numImages; i++)
	{
		int noNeighbors = matches[i].size();
		for(int j=0; j<noNeighbors; j++)
		{
			int numMatches = matches[i][j].size();
			int numValidMatches = 0;
			for (int k = 0; k < numMatches; k++) 
			{
				if( matches[i][j][k].validMatch)
				{
					numValidMatches += 1;
				}
			}
			//if (numValidMatches >= 15) 
			{
				// Write the pair 
				fprintf(f, "%d %d\n", i, i+j+1);
				
				// Write the number of matches 
				fprintf(f, "%d\n", numValidMatches);
				
				for (int k = 0; k < numMatches; k++) 
				{
					if( matches[i][j][k].validMatch)
					{
						fprintf(f, "%d %d\n",  matches[i][j][k].m_idx1,  matches[i][j][k].m_idx2);
					}
				}
			}
		}
	}
	fclose(f);
	
}

void KeypointMatcher::readMatches(const char* filename, vector<vector<vector<KeypointMatch> > > &matches)
{
	ifstream fin(filename, ios::in);
	int img1, img2;
	int noMatches;
	while(!fin.eof())
	{
		vector<KeypointMatch> keyMatches;
		fin >> img1 >> img2;
		if(img1 == 23 && img2==24)
			int debug = 1;
		
		fin >> noMatches;
		for(int m=0; m<noMatches; m++)
		{
			KeypointMatch km;
			fin >> km.m_idx1 >> km.m_idx2;
			km.validMatch = true;
			keyMatches.push_back(km);
		}
		matches[img1][img2-img1-1] = (keyMatches);
		printf("img %d to %d has %d matches\n", img1, img2, keyMatches.size());
	}
	fin.close();
}

// Compute a set of tracks that explain the matches
void KeypointMatcher::computeTracks(vector<vector<keypt_t> > &keyInfo, vector<vector<vector<KeypointMatch> > > &matches, std::vector<CorrespondenceTrack> &tracks) 
{
    unsigned int numImages = matches.size();
	
    // Clear all marks
    for (unsigned int i = 0; i < numImages; i++) 
	{
		int noFeatures = keyInfo[i].size();
		for(int j=0; j<noFeatures; j++)
		{
			keyInfo[i][j].visited = false;
		}
	}
	
	// Sort all match lists
	for (unsigned int i = 0; i < numImages; i++) 
	{
		for(int j=0; j<matches[i].size(); j++)
		{
			std::vector<KeypointMatch> &list = matches[i][j];
			sort(list.begin(), list.end(), compareFirst);
		}
    }
		
	std::vector<bool> imgMarks;
	imgMarks.resize(numImages, false);
	
	for(int i=0; i<numImages; i++)
	{
		int noNeighbors = matches[i].size();
		if(noNeighbors == 0)
			continue;
		
		int noFeatures = keyInfo[i].size();
		for(int j=0; j<noFeatures; j++)
		{
			if(keyInfo[i][j].visited)
				continue;
			
			imgMarks.clear();
			imgMarks.resize(numImages, false);
			
			std::vector<ImageKey> features;
			std::queue<ImageKey> featuresQueue;
			
			keyInfo[i][j].visited = true;
			features.push_back(ImageKey(i,j));
			featuresQueue.push(ImageKey(i,j));
			imgMarks[i] = true;
			
			int noRounds = 0;
			while(!featuresQueue.empty())
			{
				noRounds++;
				ImageKey f = featuresQueue.front();
				featuresQueue.pop();
				
				int imgIndex = f.first;
				int keyIndex = f.second;
				
				KeypointMatch tmp;
				tmp.m_idx1 = keyIndex;
				
				for(int k=0;k<matches[imgIndex].size(); k++)
				{
					int newImgIndex = imgIndex+k+1;
					
					//image visited already
					if(imgMarks[newImgIndex])
						continue;
					
					std::vector<KeypointMatch> &list = matches[imgIndex][k];
					std::pair<std::vector<KeypointMatch>::iterator, std::vector<KeypointMatch>::iterator> p;
					
                    p = equal_range(list.begin(), list.end(), tmp, compareFirst);
					
                    if(p.first == p.second)
                        continue;
					
					if((p.first)->m_idx1 != keyIndex)
						continue;
					
					int newKeyIndex = (p.first)->m_idx2;
					
					if(newKeyIndex >= keyInfo[newImgIndex].size())
						continue;
				
					//feature visited already
					if(keyInfo[newImgIndex][newKeyIndex].visited)
						continue;
					
                    // Mark and push the point
					keyInfo[newImgIndex][newKeyIndex].visited = true;
					imgMarks[newImgIndex] = true;
					features.push_back(ImageKey(newImgIndex, newKeyIndex));
					featuresQueue.push(ImageKey(newImgIndex, newKeyIndex));
				}
			}
			if (features.size() == 15) 
			{
				printf("Point with %d projections found\n", (int) features.size());
				
				CorrespondenceTrack ct;
				
				for(int k=0; k<features.size(); k++)
				{
					Correspondence c;
					keypt_t kp = keyInfo[features[k].first][features[k].second];
					keyInfo[features[k].first][features[k].second].trackIndex = tracks.size();
					c.imageIndex = features[k].first;
					c.keyIndex = features[k].second;
					c.position = Vec2f(kp.x, kp.y);
					ct.correspondences.push_back(c);
					if(features.size() == 5)
					{
						printf("In image %d:%f %f\n", c.imageIndex, c.position[0], c.position[1]);
					}
				}
				
				tracks.push_back(ct);
			}
		}
	}
}

void KeypointMatcher::setTrackMatches(std::vector<CorrespondenceTrack> &tracks, int numImages, std::vector<std::vector<std::vector<KeypointMatch> > > &trackMatches)
{
	int numTracks = tracks.size();
	
	std::vector<std::vector<std::vector<KeypointMatch> > > tmpMatches;
	trackMatches.resize(numImages);
	for(int i=0; i<numImages; i++)
	{
		trackMatches[i].resize(numImages-1-i);
	}
	
	for (int i = 0; i < numTracks; i++) 
	{
		CorrespondenceTrack ct = tracks[i];
		int numViews = ct.correspondences.size();
		
		for (int j = 0; j < numViews; j++) 
		{
			int img1 = ct.correspondences[j].imageIndex;
			
			for (int k = j+1; k < numViews; k++) 
			{
				int img2 = ct.correspondences[k].imageIndex;
				
				if(img1 < img2)
				{
					KeypointMatch kp;
					kp.m_idx1 = ct.correspondences[j].keyIndex;
					kp.m_idx2 = ct.correspondences[k].keyIndex;
					trackMatches[img1][img2-img1-1].push_back(kp);
				}
				else
				{
					KeypointMatch kp;
					kp.m_idx1 = ct.correspondences[j].keyIndex;
					kp.m_idx2 = ct.correspondences[k].keyIndex;
					trackMatches[img2][img1-img2-1].push_back(kp);
				}
			}
		}
	}
}
