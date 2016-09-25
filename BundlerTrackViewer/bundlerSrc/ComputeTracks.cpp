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

/* ComputeTracks.cpp */
/* Code for linking matches into tracks */

#include <queue>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "keys.h"

#include "BundlerApp.h"
#include "SifterUtil.h"

bool CompareFirst(const KeypointMatch &k1, const KeypointMatch &k2) {
    return (k1.m_idx1 < k2.m_idx1);
}

/* Compute a set of tracks that explain the matches */
#define LARGE_NUMBER 99999999
void BundlerApp::ComputeTracks(int new_image_start, int num_grids, int *num_grid_rows, int *num_grid_cols) 
{
    unsigned int num_images = GetNumImages();

    /* Clear all marks for new images */
    for (unsigned int i = 0; i < num_images; i++) 
	{
        // If this image has no neighbors, don't worry about its keys
        int num_nbrs = (int) m_matches.GetNumNeighbors(i);

        if (num_nbrs == 0)
            continue;

		int num_features = m_image_data[i].GetNumKeys();
        m_image_data[i].m_key_flags.resize(num_features);

	}

    int pt_idx = 0;

    std::vector<TrackData> tracks;
	std::vector<std::vector<TrackData> > gridTracks;
	std::vector<TrackData> nonGridTracks;
	
	std::vector<std::vector<bool> > gridMask;
	
	int totalGridPts = 0;
	
	if(num_grids > 0)
	{
		gridTracks.resize(num_grids);
		gridMask.resize(num_grids);
	}
	
	for(int i=0; i<num_grids; i++)
	{
		gridTracks[i].resize(num_grid_rows[i] * num_grid_cols[i]);
		gridMask[i].resize(num_grid_rows[i] * num_grid_cols[i], false);
		totalGridPts += num_grid_rows[i] * num_grid_cols[i];
	}

    // Sort all match lists
	for (unsigned int i = 0; i < num_images; i++) 
	{
        MatchAdjList::iterator iter;
        for (iter = m_matches.Begin(i); iter != m_matches.End(i); iter++) 
		{
            std::vector<KeypointMatch> &list = iter->m_match_list;
            sort(list.begin(), list.end(), CompareFirst);
        }
    }

    bool *img_marked = new bool[num_images];
    memset(img_marked, 0, num_images * sizeof(bool));

    std::vector<int> touched;
    touched.reserve(num_images);

    for (unsigned int i = 0; i < num_images; i++) 
	{
		int num_features = m_image_data[i].GetNumKeys();

        // If this image has no neighbors, skip it
        int num_nbrs = (int) m_matches.GetNumNeighbors(i);

        if (num_nbrs == 0)
            continue;

		for (int j = 0; j < num_features; j++) 
		{
			ImageKeyVector features;
			std::queue<ImageKey> features_queue;

			// Check if this feature was visited
			if (m_image_data[i].m_key_flags[j])
                continue; // already visited this feature

            // Reset flags
            int num_touched = touched.size();
            for (int k = 0; k < num_touched; k++)
                img_marked[touched[k]] = false;
            touched.clear();

			// Do a breadth first search given this feature
			m_image_data[i].m_key_flags[j] = true;

			features.push_back(ImageKey(i, j));
			features_queue.push(ImageKey(i, j));
			
			
            img_marked[i] = true;
            touched.push_back(i);

			int num_rounds = 0;
			while (!features_queue.empty()) 
			{
				num_rounds++;

				ImageKey feature = features_queue.front();
				features_queue.pop();
		
				int img1 = feature.first;
				int f1 = feature.second;
                KeypointMatch dummy;
                dummy.m_idx1 = f1;

				int start_idx;
				// Limit new images to point only to other new images
				if (img1 >= new_image_start) 
				{
					start_idx = new_image_start;
				} 
				else 
				{
					start_idx = 0;
				}

				MatchAdjList &nbrs = m_matches.GetNeighbors(img1);

                MatchAdjList::iterator iter;
                for (iter = nbrs.begin(); iter != nbrs.end(); iter++) 
				{
                    unsigned int k = iter->m_index; // *iter; // nbrs[nbr];

	
                    if (img_marked[k])
                        continue;

					MatchIndex base = GetMatchIndex(img1, k);

                    std::vector<KeypointMatch> &list = m_matches.GetMatchList(base); // m_match_lists[base];

                    // Do a binary search for the feature
                    std::pair<std::vector<KeypointMatch>::iterator, 
                              std::vector<KeypointMatch>::iterator> p;

                    p = equal_range(list.begin(), list.end(), dummy, CompareFirst);

                    if (p.first == p.second)
                        continue;  // not found

                    assert((p.first)->m_idx1 == f1);
                    int idx2 = (p.first)->m_idx2;
			    
                    // Check if we visited this point already
                    assert(idx2 < m_image_data[k].GetNumKeys());

                    if (m_image_data[k].m_key_flags[idx2])
                        continue;

                    // Mark and push the point
                    // GetKey(k,idx2).m_extra = pt_idx;
                    m_image_data[k].m_key_flags[idx2] = true;
                    features.push_back(ImageKey(k, idx2));
                    features_queue.push(ImageKey(k, idx2));

                    img_marked[k] = true;
                    touched.push_back(k);
				}
			} /* while loop */

			if (features.size() >= 2) 
			{
				//printf("Point %d has %d projections\n", features[0].second, (int) features.size()); // , num_inconsistent);
				fflush(stdout);
				
				int keyId = features[0].second;
				
				if(num_grids > 0 && keyId < totalGridPts)
				{
					int prevPts = 0;
					for(int g=0; g<num_grids; g++)
					{
						if(keyId < prevPts + num_grid_rows[g]*num_grid_cols[g])
						{
							gridTracks[g][keyId - prevPts] = TrackData(features);
							gridMask[g][keyId - prevPts] = true;
							break;
						}
						prevPts += num_grid_rows[g]*num_grid_cols[g];
					}
				}
				else if(num_grids > 0)
				{
					nonGridTracks.push_back(TrackData(features));
				}
				else
				{
					tracks.push_back(TrackData(features));
				}
				pt_idx++;
			}

		} /* for loop over features */
    } /* for loop over images */

    printf("[SifterApp::ComputeTracks] Found %d points\n", pt_idx);
    fflush(stdout);
	
	if(num_grids > 0)
	{
		//place grid tracks
		for(int g=0; g<num_grids; g++)
		{
			for(int i=0; i<num_grid_rows[g]*num_grid_cols[g]; i++)
			{
				if(gridMask[g][i])
				{
					tracks.push_back(gridTracks[g][i]);
					//printf("Grid point has %d projections\n", gridTracks[g][i].m_views.size()); 
				}
			}
		}
		//place non grid tracks
		for(int i=0; i<nonGridTracks.size(); i++)
		{
			tracks.push_back(nonGridTracks[i]);
			//printf("Point has %d projections\n", nonGridTracks[i].m_views.size()); 
		}
	}

    if (pt_idx != (int) tracks.size()) 
	{
		printf("[SifterApp::ComputeTracks] Error: point count inconsistent!\n");
		fflush(stdout);
    }

    /* Clear match lists */
    printf("[SifterApp::ComputeTracks] Clearing match lists...\n");
    fflush(stdout);

    RemoveAllMatches();

    /* Create the new consistent match lists */
    printf("[SifterApp::ComputeTracks] Creating consistent match lists...\n");
    fflush(stdout);

    int num_pts = pt_idx;

    for (int i = 0; i < num_pts; i++) 
	{
		int num_features = (int) tracks[i].m_views.size();

        for (int j = 0; j < num_features; j++) {
            int img1 = tracks[i].m_views[j].first;
            int key1 = tracks[i].m_views[j].second;

            m_image_data[img1].m_visible_points.push_back(i);
            m_image_data[img1].m_visible_keys.push_back(key1);
        }        
    }

#ifdef SBK_OUTPUT
    printf("%% Number of tracks\n");
    printf("%d\n", num_pts);

    printf("%% Number of frames\n");
    printf("%d\n", num_images);

    for (int i = 0; i < num_pts; i++) {
	int num_features = (int) tracks[i].size();

	printf("%% Length of track %d\n", i);
	printf("%d\n", num_features);
	printf("%% Track# Frame# X Y\n");

	for (int j = 0; j < num_features; j++) {
	    int img = tracks[i].m_views[j].first;
	    int f = tracks[i].m_views[j].second;
	    int w = m_image_data[img].GetWidth();
	    int h = m_image_data[img].GetHeight();

	    double x = GetKey(img,f)->m_x + 0.5 * w;
	    double y = GetKey(img,f)->m_y + 0.5 * h;

	    printf("%d %d %0.5f %0.5f\n", i, img, x, y);
	}
    }
#endif

    /* Save the tracks */
    m_track_data = tracks;

    // SetMatchesFromTracks();

    printf("[SifterApp::ComputeTracks] Done!\n");
    fflush(stdout);
}
#undef LARGE_NUMBER
