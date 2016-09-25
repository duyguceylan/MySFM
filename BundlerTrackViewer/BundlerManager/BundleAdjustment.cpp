/*
 *  BundleAdjustment.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#define MATCH_THRESHOLD 32
#define MIN_SCORE 1.0e-1
#define MIN_MATCHES 80
#define SCORE_THRESHOLD 2.0
#define INITIAL_FOCAL_LENGTH 1944.887
#define INITIAL_DEPTH 3.0

#include "BundleAdjustment.h"

BundleAdjustment::BundleAdjustment()
{
	projectionEstimationThreshold = 4.0;
}

//Pick a good initial pair of cameras to bootstrap the bundle adjustment
void BundleAdjustment::pickInitialPair(int &imgIndex1, int &imgIndex2, bool useInitFocalOnly, std::vector<std::vector<std::vector<KeypointMatch> > > &trackMatches,
									   std::vector<std::vector<Transformation> > &transformations)
{
    // Compute the match matrix
    int numImages = trackMatches.size();
	
	int maxMatches = 0;
	double maxScore = 0.0;
	double maxScoreRes = 0.0;
	
	int imgResIndex1 = -1;
	int imgResIndex2 = -1;
	imgIndex1 = -1;
	imgIndex2 = -1;
	
	int maxPts = 0;
	for(int i=0; i<numImages; i++)
	{
		if(useInitFocalOnly /*&& TODO: has focal length*/)
			continue;
		for(int j=i+1; j<numImages; j++)
		{
			if(useInitFocalOnly /*&& TODO: has focal length*/ )
				continue;
			
			int numMatches = trackMatches[i][j-i-1].size();
			maxPts += numMatches;
			
			if(numMatches < MATCH_THRESHOLD)
				continue;
			
			double score = 0.0;
			double ratio = transformations[i][j].inlierRatio;
			
			if (ratio == 0.0) 
			{
				score = MIN_SCORE;
			} 
			else 
			{
				score = 1.0 / ratio;
			}
			
			if (numMatches > maxMatches && score > SCORE_THRESHOLD) 
			{
                maxMatches = numMatches;
                maxScore = score;
				
                imgIndex1 = i;
                imgIndex2 = j;
            }
			
            if (numMatches > MIN_MATCHES && score > maxScoreRes) 
			{
                maxScoreRes = score;
                imgResIndex1 = i;
                imgResIndex2 = j;
            }
			
		}
	}
	
	if (imgIndex1 == -1 && imgIndex2 == -1) 
	{
        if (imgResIndex1 == -1 && imgResIndex2 == -1) 
		{
            printf("[BundleAdjustment] Error: no good camera pairs found!\n");
			
            if (useInitFocalOnly) 
			{
                printf("[BundleAdjustment] Trying a backup approach...\n");
                pickInitialPair(imgIndex1, imgIndex2, false, trackMatches, transformations);
            } 
			else 
			{
                printf("[BundleAdjustment] Picking first two cameras...\n");
				
                imgIndex1 = 0;
                imgIndex2 = 1;
            }
        } 
		else 
		{
            imgIndex1 = imgResIndex1;
            imgIndex2 = imgResIndex2;
        }
    }
}

void BundleAdjustment::initializeCameraParams(camera_params_t &camera)
{
	camera.R[0] = 1.0; camera.R[1] = 0.0; camera.R[2] = 0.0;
	camera.R[3] = 0.0; camera.R[4] = 1.0; camera.R[5] = 0.0;
	camera.R[6] = 0.0; camera.R[7] = 0.0; camera.R[8] = 1.0;
    
	camera.t[0] = camera.t[1] = camera.t[2] = 0.0;
    
	camera.f = 0.0;
    
	camera.k[0] = camera.k[1] = 0.0;
	
    camera.k_inv[0] = camera.k_inv[2] = camera.k_inv[3] = 0.0;
    camera.k_inv[4] = camera.k_inv[5] = 0.0;
    camera.k_inv[1] = 1.0;
	
    camera.f_scale = 1.0;
    camera.k_scale = 1.0;
	
    for (int i = 0; i < NUM_CAMERA_PARAMS; i++) 
	{
        camera.constrained[i] = 0;
        camera.constraints[i] = 0.0;
        camera.weights[i] = 0.0;
    }
	
	//TODO
    //if (data.m_known_intrinsics) 
	//{
    //    camera.known_intrinsics = 1;
    //    memcpy(camera.K_known, data.m_K, 9 * sizeof(double));
    //    memcpy(camera.k_known, data.m_k, 5 * sizeof(double));
    //}
	//else 
	//{
    //    camera.known_intrinsics = 0;
    //}
}

void BundleAdjustment::setCameraConstraints(int imgIndex, camera_params_t &cam)
{
	//TODO
}

/* Setup the initial camera pair for bundle adjustment */
void BundleAdjustment::setupInitialCameraPair(int imageIndex1, int imageIndex2, double &initFocalLength1, double &initFocalLength2,
											  std::vector<camera_params_t> &cameras, std::vector<std::vector<keypt_t> > &keyInfo, 
											 std::vector<std::vector<std::vector<KeypointMatch> > > &trackMatches, 
											 std::vector<CorrespondenceTrack> &tracks)
{
	initializeCameraParams(cameras[imageIndex1]);
	initializeCameraParams(cameras[imageIndex2]);
	
	setCameraConstraints(imageIndex1, cameras[imageIndex1]);
	setCameraConstraints(imageIndex2, cameras[imageIndex2]);
	
	//put first camera at origin
	cameras[imageIndex1].R[0] = 1.0;  cameras[imageIndex1].R[1] = 0.0;  cameras[imageIndex1].R[2] = 0.0;
    cameras[imageIndex1].R[3] = 0.0;  cameras[imageIndex1].R[4] = 1.0;  cameras[imageIndex1].R[5] = 0.0;
    cameras[imageIndex1].R[6] = 0.0;  cameras[imageIndex1].R[7] = 0.0;  cameras[imageIndex1].R[8] = 1.0;
	
	//TODO
	//if (m_image_data[i_best].m_camera.m_constrained[0])
    //    cameras[0].t[0] = m_image_data[i_best].m_camera.m_constraints[0];
    //else
    //    cameras[0].t[0] = 0.0;
    //if (m_image_data[i_best].m_camera.m_constrained[1])
    //    cameras[0].t[1] = m_image_data[i_best].m_camera.m_constraints[1];
    //else
    //    cameras[0].t[1] = 0.0;
    //if (m_image_data[i_best].m_camera.m_constrained[2])
    //    cameras[0].t[2] = m_image_data[i_best].m_camera.m_constraints[2];
    //else
    //    cameras[0].t[2] = 0.0;
	cameras[imageIndex1].t[0] = 0.0;  cameras[imageIndex1].t[1] = 0.0;	cameras[imageIndex1].t[2] = 0.0;
	
	//TODO
	//if (m_image_data[i_best].m_has_init_focal)
    //    init_focal_length_0 = cameras[0].f = m_image_data[i_best].m_init_focal;
    //else 
    //    init_focal_length_0 = cameras[0].f = m_init_focal_length; // INITIAL_FOCAL_LENGTH;
	
    //if (m_image_data[j_best].m_has_init_focal)
    //    init_focal_length_1 = cameras[1].f = m_image_data[j_best].m_init_focal;
    //else
    //    init_focal_length_1 = cameras[1].f = m_init_focal_length; // INITIAL_FOCAL_LENGTH;
	initFocalLength1 = cameras[imageIndex1].f = INITIAL_FOCAL_LENGTH;
	initFocalLength2 = cameras[imageIndex2].f = INITIAL_FOCAL_LENGTH;
	
    //TODO
	bool solvedForExtrinsics = false;
    //if (m_factor_essential && m_image_data[i_best].m_has_init_focal && 
    //    m_image_data[j_best].m_has_init_focal && 
    //    !m_use_constraints) 
	//{
	//	
	//	/* Solve for the initial locations */
	//	if (EstimateRelativePose2(i_best, j_best, cameras[0], cameras[1])) {
	//		solved_for_extrinsics = true;
	//	}        
    //} 
	//else 
	{
		// Put second camera at origin too
        cameras[imageIndex2].R[0] = 1.0;  cameras[imageIndex2].R[1] = 0.0;  cameras[imageIndex2].R[2] = 0.0;
        cameras[imageIndex2].R[3] = 0.0;  cameras[imageIndex2].R[4] = 1.0;  cameras[imageIndex2].R[5] = 0.0;
        cameras[imageIndex2].R[6] = 0.0;  cameras[imageIndex2].R[7] = 0.0;  cameras[imageIndex2].R[8] = 1.0;
		
        //TODO
		//if (m_image_data[j_best].m_camera.m_constrained[0])
        //    cameras[1].t[0] = m_image_data[j_best].m_camera.m_constraints[0];
        //else
        //    cameras[1].t[0] = 0.0;
        //if (m_image_data[j_best].m_camera.m_constrained[1])
        //    cameras[1].t[1] = m_image_data[j_best].m_camera.m_constraints[1];
        //else
        //    cameras[1].t[1] = 0.0;
        //if (m_image_data[j_best].m_camera.m_constrained[2])
        //    cameras[1].t[2] = m_image_data[j_best].m_camera.m_constraints[2];
        //else
        //    cameras[1].t[2] = 0.0;   
		cameras[imageIndex2].t[0] = 0.0;  cameras[imageIndex2].t[1] = 0.0;	cameras[imageIndex2].t[2] = 0.0;
    }
	
	//set focal constraint
	cameras[imageIndex1].constrained[6] = true;
	cameras[imageIndex1].constraints[6] = INITIAL_FOCAL_LENGTH;
	cameras[imageIndex1].weights[6] = 0.01;
	
	cameras[imageIndex2].constrained[6] = true;
	cameras[imageIndex2].constraints[6] = INITIAL_FOCAL_LENGTH;
	cameras[imageIndex2].weights[6] = 0.01;
		
    // Set up the initial 3D points
    printf("[BundleAdjustment] Adding initial matches...\n");
	
    std::vector<KeypointMatch> &list = trackMatches[imageIndex1][imageIndex2];
	
    unsigned int numMatches = list.size();
	
    for (unsigned int i = 0; i < numMatches; i++) 
	{
        int keyIdx1 = list[i].m_idx1;
        int keyIdx2 = list[i].m_idx2;
		
        printf("  Adding match %d ==> %d\n", keyIdx1, keyIdx2);
		
        double xProj = keyInfo[imageIndex1][keyIdx1].x;
        double yProj = keyInfo[imageIndex1][keyIdx1].y;
		
        // Back project the point to a constant depth
        if (!solvedForExtrinsics) 
		{
            double xPt = (xProj / cameras[imageIndex1].f) * INITIAL_DEPTH;
            double yPt = (yProj / cameras[imageIndex1].f) * INITIAL_DEPTH;
            double zPt = INITIAL_DEPTH + cameras[imageIndex1].t[2];
			
			tracks[keyInfo[imageIndex1][keyIdx1].trackIndex].pos3D = Vec3f(xPt, yPt, zPt);
        } 
		else 
		{
            double xProj1 = keyInfo[imageIndex1][keyIdx1].x;
			double yProj1 = keyInfo[imageIndex1][keyIdx1].y;
			double xProj2 = keyInfo[imageIndex2][keyIdx2].x;
			double yProj2 = keyInfo[imageIndex2][keyIdx2].y;
			
            double error;
			
            Vec2f p(xProj1, yProj1);
            Vec2f q(xProj2, yProj2);
			
            bool inFront = true;
            double angle = 0.0;
            Vec3f pt;// = Triangulate(p, q, cameras[imageIndex1], cameras[imageIndex2], error, inFront, angle, true);
			
            printf(" tri.error[%d] = %0.3f\n", i, error);
			
            if (error > projectionEstimationThreshold) 
			{
                printf(" skipping point\n");
                continue;
            }
			else
			{
				tracks[keyInfo[imageIndex1][keyIdx1].trackIndex].pos3D = pt;
			}
        }
		
        //TODO
		// Get the color of the point
        //unsigned char r = GetKey(i_best,key_idx1).m_r;
        //unsigned char g = GetKey(i_best,key_idx1).m_g;
        //unsigned char b = GetKey(i_best,key_idx1).m_b;
        //colors[pt_count] = v3_new((double) r, (double) g, (double) b);
		
	}
}
