/*
 *  EpipolarGeometry.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "EpipolarGeometry.h"
#include "../MathUtils/MathUtils.h"
#include "../MathUtils/GeometryUtils.h"
#include "Line.h"
#include "vector.h"
#include "imagelib/fmatrix.h"
#include "5point/5point.h"
#include "Mesh3D.h"

float EpipolarGeometry::globalScale = 1.0;
std::vector<Vec3f> EpipolarGeometry::globalOuts = vector<Vec3f> ();
std::vector<Vec3f> EpipolarGeometry::globalIns = vector<Vec3f> ();

void EpipolarGeometry::fmatrixResiduals(const int *m, const int *n, const double *x, double *fvec, int *iflag) 
{
    int i;
    double sum = 0.0;
	
	Matrix3f F;
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			if(i==2 && j==2)
				break;
			
			F[i][j] = x[i*3+j];
		}
	}
	
    F[2][2] = EpipolarGeometry::globalScale;
	
	Matrix3f F2;
	MathUtils::closestRank2Matrix(F, F2);
	
	int globalNumMatches = globalOuts.size();
	
    if (globalNumMatches != (*m)) 
	{
		printf("Error: number of matches don't match!\n");
		return;
    }
	
    for (i = 0; i < globalNumMatches; i++) 
	{
		fvec[i] = sqrt(fmatrixComputeResidual(F2, globalOuts[i], globalIns[i]));
		if (*iflag == 0) 
		{
			sum += fvec[i];
		}
    }
}

float EpipolarGeometry::fmatrixComputeResidual(Matrix3f F, Vec3f l, Vec3f r)
{
	Vec3f Fl = F.transpose() * l;
	Vec3f Fr = F * r;
	
	float pt = dot(l, Fr);
	
	float resid = (1.0 / ((Fl[0] * Fl[0] + Fl[1] * Fl[1]) + (Fr[0] * Fr[0] + Fr[1] * Fr[1]))) * (pt * pt);
	return resid;
}

int EpipolarGeometry::estimateFmatrixLinear(std::vector<Vec3f> &rightPoints, std::vector<Vec3f> &leftPoints, Matrix3f &Fout) 
{
	int numPts = rightPoints.size();
	if(numPts < 8)
	{
		printf("[estimateFmatrixLinear] Insufficient correspondences (need at least 8, given only %d)\n", numPts);
		return 0;
	}
	
	//First, compute the centroid of both sets of points
	Vec3f rightCentroid(0.0, 0.0, 0.0);
	Vec3f leftCentroid(0.0, 0.0, 0.0);
	for(int i=0; i<numPts; i++)
	{
		rightCentroid += rightPoints[i];
		leftCentroid += leftPoints[i];
	}
	
	rightCentroid /= numPts;
	leftCentroid /= numPts;
	
	// Compute the average distance from each point to the centroid
    float rDist = 0;
	float lDist = 0;
    
    for (int i = 0; i < numPts; i++) 
	{
		rDist += (rightCentroid - rightPoints[i]).length();
		lDist += (leftCentroid - leftPoints[i]).length();
    }
	
    rDist /= numPts;
    lDist /= numPts;

	//scale the avg distance
	float rScale = sqrt(2.0) / rDist;
	float lScale = sqrt(2.0) / lDist;
	
	//normalize the points
	Matrix3f H, Hp;
	H.setIdentity();
	Hp.setIdentity();
	
	H[0][0] = lScale;
	H[2][0] = -lScale * leftCentroid[0];
	H[1][1] = lScale;
	H[2][1] = -lScale*leftCentroid[1];
	
	Hp[0][0] = rScale;
	Hp[2][0] = -rScale * rightCentroid[0];
	Hp[1][1] = rScale;
	Hp[2][1] = -rScale*rightCentroid[1];
	
	vector<Vec3f> newRightPoints, newLeftPoints;
	
	for(int i=0; i<numPts; i++)
	{
		newRightPoints.push_back(Hp*rightPoints[i]);
		newLeftPoints.push_back(H*leftPoints[i]);
		//printf("l:%f %f %f, r:%f %f %f\n", newLeftPoints[i][0], newLeftPoints[i][1], newLeftPoints[i][2], newRightPoints[i][0], newRightPoints[i][1], newRightPoints[i][2]);
	}
	
	
	//fill in matrix A
	float *A = new float[numPts*8];
	float *b = new float[numPts];
	for(int i=0; i<numPts; i++)
	{
		A[i] = newLeftPoints[i][0]*newRightPoints[i][0];
		A[numPts + i] = newLeftPoints[i][0]*newRightPoints[i][1];
		A[numPts*2 + i] = newLeftPoints[i][0];
		A[numPts*3 + i] = newLeftPoints[i][1]*newRightPoints[i][0];
		A[numPts*4 + i] = newLeftPoints[i][1]*newRightPoints[i][1];
		A[numPts*5 + i] = newLeftPoints[i][1];
		A[numPts*6 + i] = newRightPoints[i][0];
		A[numPts*7 + i] = newRightPoints[i][1];
		b[i] = -1.0;
		//A[numPts*8 + i] = 1.0;
		//printf("%f %f %f %f %f %f %f %f %f\n", A[i], A[numPts + i], A[numPts*2 + i], A[numPts*3 + i], A[numPts*4 + i], A[numPts*5 + i], A[numPts*6 + i],
		//	   A[numPts*7 + i], A[numPts*8 + i]);
	}
	
	int noRows = numPts;
	int noCols = 1;//9;
	int success;
	
	success= MathUtils::solveLinearSystem(A, b, &noRows, &noCols);
	
	double *u = new double[noRows*noRows];
	double *s = new double[noRows*noCols];
	double *vt = new double[noCols*noCols];
	memset(u, 0, sizeof(double)*noRows*noRows);
	memset(s, 0, sizeof(double)*noRows*noCols);
	memset(vt, 0, sizeof(double)*noCols*noCols);
	//success = MathUtils::computeSVD(A, &noRows, &noCols, u, s, vt);
	
	if(success != 0)
	{
		printf("[estimateFmatrixLinear] could not compute SVD.\n");
		return 0;
	}
    
	//find F as the last column (last row) of v(vt) in row major
	//needs to be converted to column major for second call to svd
	Matrix3f F;
	F[0][0] = b[0]; F[0][1] = b[1]; F[0][2] = b[2];
	F[1][0] = b[3]; F[1][1] = b[4]; F[1][2] = b[5];
	F[2][0] = b[6]; F[2][1] = b[7]; F[2][2] = 1.0;
	
	double *lastRowVt = new double[9];
	for(int i=0; i<9; i++)
	{
		lastRowVt[i] = vt[i*9+8];
	}
	
	//for(int i=0; i<3; i++)
	//{
	//	for(int j=0; j<3; j++)
	//	{
	//		F[i][j] = lastRowVt[i*3+j]; 
	//	}
	//}
	
	// Un-normalize
	F = Hp.transpose() * F * H;
	
	//constrain F to have determinant 0
	MathUtils::closestRank2Matrix(F, Fout);
	
    //delete [] u;
	//delete [] s;
	//delete [] vt;
	//delete [] lastRowVt;
	return 1;
}

void EpipolarGeometry::refineFmatrixNonlinearMatches(int numPts, std::vector<Vec3f> &rightPoints, std::vector<Vec3f> &leftPoints, Matrix3f &F, Matrix3f &FOut)
{
    double Ftmp[9];
	
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			Ftmp[i*3+j] = F[i][j];
		}
	}
	
    globalScale = F[2][2];
	
	globalIns.clear();
	globalIns.assign(leftPoints.begin(), leftPoints.end());
	
	globalOuts.clear();
	globalOuts.assign(rightPoints.begin(), rightPoints.end());
	
	MathUtils::lmdifDriver2((EpipolarGeometry::fmatrixResiduals), numPts, 8, Ftmp, 1.0e-12);
	
    Ftmp[8] = EpipolarGeometry::globalScale;
	//printf("FOut:%f %f %f %f %f %f %f %f %f\n", Ftmp[0], Ftmp[1], Ftmp[2], Ftmp[3], Ftmp[4], Ftmp[5], Ftmp[6], Ftmp[7], Ftmp[8]);
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			FOut[i][j] = Ftmp[i*3+j];
		}
	}
	
	Matrix3f F2;
	MathUtils::closestRank2Matrix(FOut, F2);
	FOut = F2;
}

// Use RANSAC to estimate an F-matrix
int EpipolarGeometry::estimateFmatrixRansacMatches(int numPts, vector<Vec3f> &aPts, vector<Vec3f> &bPts, int numTrials, double threshold, 
												   double successRatio, Matrix3f &F, vector<int> &inliers) 
{
	if (numPts < 8) 
	{
		printf("[estimateFmatrixRansacMatches] Could not find 8 good correspondences, F-matrix estimation failed\n");
		return 0;
    }

	int success;
	int numInliersMax = 0;
	float totalDist = 0.0;
	Matrix3f FBest;
	
	for(int i=0; i<numTrials; i++)
	{
		int idxs[8];
		int round = 0;
		
		// Sample 8 random correspondences
		for (int j = 0; j < 8; j++) 
		{
			int reselect = 0;
			
            if (round == 1000)
                return 0;
			
			int idx = rand() % numPts;
			
			/* Make sure we didn't sample this index yet */
			for (int k = 0; k < j; k++) 
			{
				if (idx == idxs[k] || aPts[idx] == aPts[idxs[k]] || bPts[idx] == bPts[idxs[k]])
				{
					reselect = 1;
					break;
				}
			}
			
			if (reselect) 
			{
                round++;
				j--;
				continue;
			}
			
			idxs[j] = idx;
		}
		
		std::vector<Vec3f> leftPoints;
		std::vector<Vec3f> rightPoints;
		for(int j=0; j<8; j++)
		{
			leftPoints.push_back(bPts[idxs[j]]);
			rightPoints.push_back(aPts[idxs[j]]);
		}
		
		Matrix3f Ftmp;
		estimateFmatrixLinear(rightPoints, leftPoints, Ftmp);
		
		int numInliers = 0;
		if(success)
		{
			if(!isnan(Ftmp[0][0]) && !isnan(Ftmp[0][1]) && !isnan(Ftmp[0][2]) && !isnan(Ftmp[1][0]) && !isnan(Ftmp[1][1]) && 
			   !isnan(Ftmp[1][2]) && !isnan(Ftmp[2][0]) && !isnan(Ftmp[2][1]) && !isnan(Ftmp[2][2]))
			{
				// Compute residuals
				totalDist = 0.0;
				vector<int> tmpInliers;
				
				for (int j = 0; j < numPts; j++) 
				{
					float resid = fmatrixComputeResidual(Ftmp, aPts[j], bPts[j]);
					totalDist += resid;
					if (resid < threshold)
					{
						numInliers++;
						tmpInliers.push_back(j);
					}
				}
				
				if (numInliers > numInliersMax) 
				{
					numInliersMax = numInliers;
					FBest = Ftmp;
					inliers.clear();
					for(int j=0; j<tmpInliers.size(); j++)
					{
						inliers.push_back(tmpInliers[j]);
					}
				}
			}
		}
		
		if ((double) numInliers / numPts > successRatio)
		{
			printf("s:%d %d %d %d %d %d %d %d\n", idxs[0], idxs[1], idxs[2], idxs[3], idxs[4], idxs[5], idxs[6], idxs[7]);
			break;
		}
	}
	
	F = FBest;
	
	
	return numInliersMax;
    
}


std::vector<int> EpipolarGeometry::estimateFMatrix(std::vector<keypt_t> &k1, std::vector<keypt_t> &k2, std::vector<KeypointMatch> matches, 
												   int numTrials, double threshold, Matrix3f &F)
{
	int numPts = (int) matches.size();
	
    // num_pts should be greater than a threshold
    if (numPts < 12) 
	{
        std::vector<int> inliers;
        return inliers;
    }
	
	std::vector<Vec3f> k1Pts;
	std::vector<Vec3f> k2Pts;
	std::vector<Vec3f> k1PtsIn;
	std::vector<Vec3f> k2PtsIn;
	
    for (int i = 0; i < numPts; i++) 
	{
        int idx1 = matches[i].m_idx1;
        int idx2 = matches[i].m_idx2;
		
		k1Pts.push_back(Vec3f(k1[idx1].x, k1[idx1].y, 1.0));
		k2Pts.push_back(Vec3f(k2[idx2].x, k2[idx2].y, 1.0));
    }
	
    std::vector<int> inliers;
	
	estimateFmatrixRansacMatches(numPts, k2Pts, k1Pts, numTrials, threshold, 0.95, F, inliers);
	
	inliers.clear();
	
	double totalDist = 0;
    for (int i = 0; i < numPts; i++) 
	{
        double dist = fmatrixComputeResidual(F, k2Pts[i], k1Pts[i]);
        totalDist += dist;
		if (dist < threshold) 
		{
            inliers.push_back(i);
        }
    }
	
    //Re-estimate using inliers
    int numInliers = (int) inliers.size();
	
    Matrix3f F0;
	
	printf("inliers before:%d\n", numInliers);
	printf("totalDist before:%f\n", totalDist);
	
    //Refine using NLLS
	for (int i = 0; i < numInliers; i++)
	{
		k1PtsIn.push_back(k1Pts[inliers[i]]);
		k2PtsIn.push_back(k2Pts[inliers[i]]);
	}
		
	F0 = F;
	refineFmatrixNonlinearMatches(numInliers, k1PtsIn, k2PtsIn, F0, F);

	inliers.clear();
	totalDist = 0;
	
	for (int i = 0; i < numPts; i++) 
	{
        double dist = fmatrixComputeResidual(F, k2Pts[i], k1Pts[i]);
        totalDist += dist;
		if (dist < threshold) 
		{
            inliers.push_back(i);
        }
    }
	
    numInliers = (int) inliers.size();
	printf("inliers after:%d\n", numInliers);
	printf("totalDist after:%f\n", totalDist);
	
	k1Pts.clear();
	k2Pts.clear();
	k1PtsIn.clear();
	k2PtsIn.clear();
	
    return inliers;
}

std::vector<int> EpipolarGeometry::estimateFMatrix(int numPts, std::vector<Vec3f> &k1Pts, std::vector<Vec3f> &k2Pts, 
												   int numTrials, double threshold, Matrix3f &F)
{
	// num_pts should be greater than a threshold
    if (numPts < 16) 
	{
        std::vector<int> inliers;
        return inliers;
    }
	
	std::vector<Vec3f> k1PtsIn;
	std::vector<Vec3f> k2PtsIn;
	
    std::vector<int> inliers;
	
	estimateFmatrixRansacMatches(numPts, k2Pts, k1Pts, numTrials, threshold, 0.95, F, inliers);
	
	double totalDist = 0;
    for (int i = 0; i < numPts; i++) 
	{
        double dist = fmatrixComputeResidual(F, k2Pts[i], k1Pts[i]);
        totalDist += dist;
		if (dist < threshold) 
		{
            inliers.push_back(i);
        }
    }
	
    //Re-estimate using inliers
    int numInliers = (int) inliers.size();
	
    Matrix3f F0;
	
	printf("inliers before:%d\n", numInliers);
	printf("totalDist before:%f\n", totalDist);
	
    //Refine using NLLS
	for (int i = 0; i < numInliers; i++)
	{
		k1PtsIn.push_back(k1Pts[inliers[i]]);
		k2PtsIn.push_back(k2Pts[inliers[i]]);
	}
		
	F0 = F;
	refineFmatrixNonlinearMatches(numInliers, k2PtsIn, k1PtsIn, F0, F);

	inliers.clear();
	totalDist = 0;
	
	for (int i = 0; i < numPts; i++) 
	{
        double dist = fmatrixComputeResidual(F, k2Pts[i], k1Pts[i]);
        totalDist += dist;
		if (dist < threshold) 
		{
            inliers.push_back(i);
        }
    }
	
    numInliers = (int) inliers.size();
	printf("inliers after:%d\n", numInliers);
	printf("totalDist after:%f\n", totalDist);
	
	k1PtsIn.clear();
	k2PtsIn.clear();
	
	printf("F:\n");
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			printf("%f ", F[j][i]);
		}
		printf("\n");
	}
	
    return inliers;
}

bool EpipolarGeometry::areFMatricesEqual(Matrix3f F1, Matrix3f F2, int width, int height, vector<Vec2f> &pts, float threshold)
{
	int noPts = pts.size();
	float avgDistance = 0;
	int count = 0;
	for(int i=0; i<noPts; i++)
	{
		Vec3f point = expand2To3(pts[i]);
		Vec3f lineEq1 = F1*point;
		Vec3f lineEq2 = F2*point;
		Vec2f minCoord, maxCoord;
		if(GeometryUtils::clipLineWrt2DRectangle(lineEq1, 0, width, 0, height, minCoord, maxCoord))
		{
			Line line1(minCoord, maxCoord);
			if(GeometryUtils::clipLineWrt2DRectangle(lineEq2, 0, width, 0, height, minCoord, maxCoord))
			{
				Line line2(minCoord, maxCoord);
				Vec3f point1 = expand2To3(line1.generateRandomPointOnLine2D());
				Vec3f point2 = expand2To3(line2.generateRandomPointOnLine2D());
		
				float dist1 = line2.findPerpendicularDistanceTo2DLine(shrink3To2(point1));
				Vec3f lineEq3 = F2.transpose() * point1;
				Vec2f normal;
				float l = sqrt(lineEq3[0]*lineEq3[0] + lineEq3[1]*lineEq3[1]);
				normal[0] = lineEq3[0] / l; normal[1] = lineEq3[1] / l;
				float d = lineEq3[2] / l;
				float dist1Prime = dot(normal, shrink3To2(point)) + d;
				
				Vec3f lineEq = F1.transpose() * point1;
				l = sqrt(lineEq[0]*lineEq[0] + lineEq[1]*lineEq[1]);
				normal[0] = lineEq[0] / l; normal[1] = lineEq[1] / l;
				d = lineEq[2] / l;
				float dist = dot(normal, shrink3To2(point)) + d;
				
				float dist2 = line1.findPerpendicularDistanceTo2DLine(shrink3To2(point2));
				Vec3f lineEq4 = F1.transpose() * point2;
				l = sqrt(lineEq4[0]*lineEq4[0] + lineEq4[1]*lineEq4[1]);
				normal[0] = lineEq4[0] / l; normal[1] = lineEq4[1] / l;
				d = lineEq4[2] / l;
				float dist2Prime = dot(normal, shrink3To2(point)) + d;
		
				avgDistance += (dist1-dist1Prime)*(dist1-dist1Prime) + (dist2-dist2Prime)*(dist2-dist2Prime);
				count += 2;
			}
		}
	}
	
	avgDistance /= (float)count;
	if(avgDistance <= threshold)
		return true;
	else
		return false;
}

bool EpipolarGeometry::estimateHomography(int numPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts, Matrix3f &H)
{
	double *A = new double[numPts*8*2];
	double *b = new double[numPts*2];
	
	for(int i=0; i<numPts; i++)
	{
		A[i*2+0] = aPts[i][0];
		A[i*2+1] = 0.0;
		A[numPts*2+i*2+0] = aPts[i][1];
		A[numPts*2+i*2+1] = 0.0;
		A[numPts*4+i*2+0] = 1.0;
		A[numPts*4+i*2+1] = 0.0;
		A[numPts*6+i*2+0] = 0.0;
		A[numPts*6+i*2+1] = aPts[i][0];
		A[numPts*8+i*2+0] = 0.0;
		A[numPts*8+i*2+1] = aPts[i][1];
		A[numPts*10+i*2+0] = 0.0;
		A[numPts*10+i*2+1] = 1.0;
		A[numPts*12+i*2+0] = -aPts[i][0]*bPts[i][0];
		A[numPts*12+i*2+1] = -aPts[i][0]*bPts[i][1];
		A[numPts*14+i*2+0] = -aPts[i][1]*bPts[i][0];
		A[numPts*14+i*2+1] = -aPts[i][1]*bPts[i][1];
		b[i*2+0] = bPts[i][0];
		b[i*2+1] = bPts[i][1];
	}
	
	int noRowsA = numPts*2;
	int noColsA = 8;
	int noColsB = 1;
	int success = MathUtils::solveLeastSquares(A, b, &noRowsA, &noColsA, &noColsB);
	if(success == 0)
	{
		H[0][0] = b[0];
		H[1][0] = b[1];
		H[2][0] = b[2];
		H[0][1] = b[3];
		H[1][1] = b[4];
		H[2][1] = b[5];
		H[0][2] = b[6];
		H[1][2] = b[7];
		H[2][2] = 1.0;
	}
	else
	{
		printf("Homography estimation failed!\n");
	}
	
	for(int i=0; i<numPts; i++)
	{
		Vec3f pos = H * expand2To3(aPts[i]);
		pos[0] = pos[0]/pos[2];
		pos[1] = pos[1]/pos[2];
		float dist = (bPts[i] - shrink3To2(pos)).length();
		printf("dist:%f\n", dist);
	}
}

void EpipolarGeometry::estimateRectifiedTransformed(int numPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts,  float &scale, Vec2f &translation) 
{
	int success;
	
	double *A = new double[numPts*6];
	double *b = new double[numPts*2];
	int noRowsA = numPts*2;
	int noColsA = 3;
	int noColsB = 1;
	
	for(int j=0; j<numPts; j++)
	{
		A[j*2] = aPts[j][0];
		A[j*2+1] = aPts[j][1];
		A[numPts*2+j*2] = 1.0;
		A[numPts*2+j*2+1] = 0.0;
		A[numPts*4+j*2] = 0.0;
		A[numPts*4+j*2+1] = 1.0;
			
		b[j*2] = bPts[j][0];
		b[j*2+1] = bPts[j][1];
	}
		
	success = MathUtils::solveLeastSquares(A, b, &noRowsA, &noColsA, &noColsB);
	if(success == 0)
	{
		if(!isnan(b[0]) && !isnan(b[1]) && !isnan(b[2]))
		{
			scale = b[0];
			translation = Vec2f(b[1], b[2]);
			
		}
		else
		{
			printf("scale and translation estimation failed.\n");
		}
	}
	else
	{
		printf("scale and translation estimation failed.\n");
	}	
}


void EpipolarGeometry::estimateRectifiedTransformedRansacMatches(int numPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts, int numTrials, double threshold, 
																 vector<float> &scale, vector<Vec2f> &translation, vector<vector<int> > &inliers) 
{
	if (numPts < 3) 
	{
		printf("[estimateRectifiedTransformedRansacMatches] Could not find 3 good correspondences, estimation failed\n");
		return ;
    }
	
	int success;
	
	double *A = new double[4*3];
	double *b = new double[4];
	int noRowsA = 4;
	int noColsA = 3;
	int noColsB = 1;
	
	int noPtsUsed = 2;
	
	for(int i=0; i<numTrials; i++)
	{
		int *idxs = new int[noPtsUsed];
		int round = 0;
		
		// Sample 2 random correspondences
		for (int j = 0; j < noPtsUsed; j++) 
		{
			int reselect = 0;
			
            if (round == 1000)
                return ;
			
			int idx = rand() % numPts;
			
			/* Make sure we didn't sample this index yet */
			for (int k = 0; k < j; k++) 
			{
				if (idx == idxs[k] || aPts[idx] == aPts[idxs[k]] || bPts[idx] == bPts[idxs[k]])
				{
					reselect = 1;
					break;
				}
			}
			
			if (reselect) 
			{
                round++;
				j--;
				continue;
			}
			
			idxs[j] = idx;
		}
		
		for(int j=0; j<noPtsUsed; j++)
		{
			A[j*2] = aPts[idxs[j]][0];
			A[j*2+1] = aPts[idxs[j]][1];
			A[noPtsUsed*2+j*2] = 1.0;
			A[noPtsUsed*2+j*2+1] = 0.0;
			A[noPtsUsed*4+j*2] = 0.0;
			A[noPtsUsed*4+j*2+1] = 1.0;
			
			b[j*2] = bPts[idxs[j]][0];
			b[j*2+1] = bPts[idxs[j]][1];
		}
		
		success = MathUtils::solveLeastSquares(A, b, &noRowsA, &noColsA, &noColsB);
		int numInliers = 0;
		if(success == 0)
		{
			if(!isnan(b[0]) && !isnan(b[1]) && !isnan(b[2]))
			{
				scale.push_back(b[0]);
				translation.push_back(Vec2f(b[1], b[2]));
				
				// Compute residuals
				Matrix3f M;
				M.setIdentity();
				M[0][0] = b[0];
				M[1][1] = b[0];
				M[2][0] = b[1];
				M[2][1] = b[2];
				vector<int> indices;
				for (int j = 0; j < numPts; j++) 
				{
					float resid = (bPts[j] - shrink3To2(M*expand2To3(aPts[j]))).length();
					//printf("resid:%f\n", resid);
					if (resid < threshold)
						indices.push_back(j);
				}
				inliers.push_back(indices);
				//printf("inliers:%d\n", indices.size());
			}
		}
	}
	
}

int EpipolarGeometry::estimateRectifiedTransformedRansacMatches(int numPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts, int numTrials, double threshold, 
																double successRatio, float &scale, Vec2f &translation, vector<int> &inliers) 
{
	if (numPts < 3) 
	{
		printf("[estimateRectifiedTransformedRansacMatches] Could not find 3 good correspondences, estimation failed\n");
		return 0;
    }
	
	int success;
	int numInliersMax = 0;
	float totalDist = 0.0;
	float scaleBest;
	Vec2f translationBest;
	
	double *A = new double[4*3];
	double *b = new double[4];
	int noRowsA = 4;
	int noColsA = 3;
	int noColsB = 1;
	
	int noPtsUsed = 2;
	
	for(int i=0; i<numTrials; i++)
	{
		int *idxs = new int[noPtsUsed];
		int round = 0;
		
		// Sample 2 random correspondences
		for (int j = 0; j < noPtsUsed; j++) 
		{
			int reselect = 0;
			
            if (round == 1000)
                return 0;
			
			int idx = rand() % numPts;
			
			/* Make sure we didn't sample this index yet */
			for (int k = 0; k < j; k++) 
			{
				if (idx == idxs[k] || aPts[idx] == aPts[idxs[k]] || bPts[idx] == bPts[idxs[k]])
				{
					reselect = 1;
					break;
				}
			}
			
			if (reselect) 
			{
                round++;
				j--;
				continue;
			}
			
			idxs[j] = idx;
		}
		
		for(int j=0; j<noPtsUsed; j++)
		{
			A[j*2] = aPts[idxs[j]][0];
			A[j*2+1] = aPts[idxs[j]][1];
			A[noPtsUsed*2+j*2] = 1.0;
			A[noPtsUsed*2+j*2+1] = 0.0;
			A[noPtsUsed*4+j*2] = 0.0;
			A[noPtsUsed*4+j*2+1] = 1.0;
			
			b[j*2] = bPts[idxs[j]][0];
			b[j*2+1] = bPts[idxs[j]][1];
		}
		
		success = MathUtils::solveLeastSquares(A, b, &noRowsA, &noColsA, &noColsB);
		int numInliers = 0;
		if(success == 0)
		{
			vector<int> tmpInliers;
			
			if(!isnan(b[0]) && !isnan(b[1]) && !isnan(b[2]))
			{
				// Compute residuals
				Matrix3f M;
				M.setIdentity();
				M[0][0] = b[0];
				M[1][1] = b[0];
				M[2][0] = b[1];
				M[2][1] = b[2];
				totalDist = 0.0;
				for (int j = 0; j < numPts; j++) 
				{
					float resid = (bPts[j] - shrink3To2(M*expand2To3(aPts[j]))).length();
					//printf("resid:%f\n", resid);
					totalDist += resid;
					if (resid < threshold)
					{
						tmpInliers.push_back(j);
						numInliers++;
					}
				}
				
				if (numInliers > numInliersMax) 
				{
					//printf("inliers:\n");
					//for(int j=0; j<noPtsUsed; j++)
					//{
					//	printf("%d ", );
					//}
					numInliersMax = numInliers;
					inliers.clear();
					for(int j=0; j<tmpInliers.size(); j++)
					{
						inliers.push_back(tmpInliers[j]);
					}
					
					scaleBest = b[0];
					translationBest[0] = b[1];
					translationBest[1] = b[2];
				}
			}
		}
	}
	
	scale = scaleBest;
	translation = translationBest;
	
	printf("%d inliers out of %d points, scale:%f, translation:%f %f\n", numInliersMax, numPts, scale, translation[0], translation[1]);
	return numInliersMax;
    
}

int EpipolarGeometry::estimatePose(int noPts, vector<Vec2f> &aPts, vector<Vec2f> &bPts, Matrix3f K1, Matrix3f K2, Matrix3f &R, Vec3f &t, vector<int> &inliers)
{
	if(noPts < 20)
		return 0;
	
	double *aVec = new double[noPts*2];
	double *bVec = new double[noPts*2];
	
	for(int i=0; i<noPts; i++)
	{
		aVec[i*2] = aPts[i][0];
		aVec[i*2+1] = aPts[i][1];
		
		bVec[i*2] = bPts[i][0];
		bVec[i*2+1] = bPts[i][1];
	}
	
	double *K1Mat = new double[9];
	double *K2Mat = new double[9];
	
	K1Mat[0] = K1[0][0]; K1Mat[1] = K1[1][0]; K1Mat[2] = K1[2][0];
	K1Mat[3] = K1[0][1]; K1Mat[4] = K1[1][1]; K1Mat[5] = K1[2][1];
	K1Mat[6] = K1[0][2]; K1Mat[7] = K1[1][2]; K1Mat[8] = K1[2][2];
	
	K2Mat[0] = K2[0][0]; K2Mat[1] = K2[1][0]; K2Mat[2] = K2[2][0];
	K2Mat[3] = K2[0][1]; K2Mat[4] = K2[1][1]; K2Mat[5] = K2[2][1];
	K2Mat[6] = K2[0][2]; K2Mat[7] = K2[1][2]; K2Mat[8] = K2[2][2];
	
	double *RMat = new double[9];
	double *tVec = new double[3];
	int *inlierIndices = new int[noPts];
	
	int noInliers = compute_pose_ransac_driver(noPts, aVec, bVec, K1Mat, K2Mat, 1.25, 512, RMat, tVec, inlierIndices);
	//printf("NoInliers:%d\n", noInliers);
	
	R[0][0] = RMat[0]; R[1][0] = RMat[1]; R[2][0] = RMat[2];
	R[0][1] = RMat[3]; R[1][1] = RMat[4]; R[2][1] = RMat[5];
	R[0][2] = RMat[6]; R[1][2] = RMat[7]; R[2][2] = RMat[8];
	
	t[0] = tVec[0]; t[1] = tVec[1]; t[2] = tVec[2];
	
	for(int i=0; i<noInliers; i++)
		inliers.push_back(inlierIndices[i]);

	return noInliers;
}

/* Given an E matrix and a point correspondence, find R and t */
bool EpipolarGeometry::findRotationAndTranslationFromEssentialMatrix(Matrix3f E, Matrix3f K, Vec2f p1, Vec2f p2, Matrix3f &R, Vec3f &t)
{	
	// Now find the SVD of E
	int noRows = 3;
	int noCols = 3;
	double *U = new double[9];
	double *S = new double[9];
	double *VT = new double[9];
	double *Er = new double[9];
	memset(U, 0.0, sizeof(double)*9);
	memset(S, 0.0, sizeof(double)*9);
	memset(VT, 0.0, sizeof(double)*9);
	
	Er[0] = E[0][0]; Er[1] = E[0][1]; Er[2] = E[0][2];
	Er[3] = E[1][0]; Er[4] = E[1][1]; Er[5] = E[1][2];
	Er[6] = E[2][0]; Er[7] = E[2][1]; Er[8] = E[2][2];
	
	MathUtils::computeSVD(Er, &noRows, &noCols, U, S, VT);
	printf("S:%f %f %f %f %f %f %f %f %f\n", S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8]);
	
	Matrix3f D, DT;
	D.setZero();
	D[0][1] = -1.0;
	D[1][0] = 1.0;
	D[2][2] = 1.0;
	DT = D.transpose();
	
	Matrix3f I;
	I.setIdentity();
    
	Vec3f t0(0.0, 0.0, 0.0);
	
    // Now find R and t
	Matrix3f UMat, VtMat;
	UMat[0][0] = U[0]; UMat[0][1] = U[1]; UMat[0][2] = U[2];
	UMat[1][0] = U[3]; UMat[1][1] = U[4]; UMat[1][2] = U[5];
	UMat[2][0] = U[6]; UMat[2][1] = U[7]; UMat[2][2] = U[8];
	
	VtMat[0][0] = VT[0]; VtMat[0][1] = VT[1]; VtMat[0][2] = VT[2];
	VtMat[1][0] = VT[3]; VtMat[1][1] = VT[4]; VtMat[1][2] = VT[5];
	VtMat[2][0] = VT[6]; VtMat[2][1] = VT[7]; VtMat[2][2] = VT[8];
	
	Vec3f tu = UMat[2];
	
	Matrix3f Ra = UMat * D * VtMat;
	Matrix3f Rb = UMat * DT * VtMat;
	
    if(Ra.getDeterminant() < 0.0)
		Ra = Ra * -1.0;
	
	if(Rb.getDeterminant() < 0.0)
		Rb = Rb * -1.0;
	
	
    // Figure out which configuration is correct using the supplied point
	float error;
    Vec3f Q = triangulate(p1, p2, K, I, t0, Ra, tu, error);
	Vec3f PQ = Ra*Q + tu;
	
	if(Q[2] < PQ[2])
	{
		R = Ra;
		t = tu;
		return true;
	}
	else if(Q[2]>0.0 && PQ[2]>0.0)
	{
		R = Ra;
		t = -tu;
		return true;
	}
	else
	{
		Q = triangulate(p1, p2, K, I, t0, Rb, tu, error);
		PQ = Rb*Q + tu;
		if(Q[2] < PQ[2])
		{
			R = Rb;
			t = tu;
			return true;
		}
		else if(Q[2]>0.0 && PQ[2]>0.0)
		{
			R = Rb;
			t = -tu;
			return true;
		}
		else
		{
			printf("Error! No case found!\n");
			return false;
		}
	}
}

Vec3f EpipolarGeometry::triangulate(Vec2f p, Vec2f q, Matrix3f K, Matrix3f R0, Vec3f t0, Matrix3f R1, Vec3f t1, float &error) 
{
    double A[12];
    double b[4];
	
	Matrix3f Kinv = K.getInverseMatrix();
	
	p = shrink3To2(Kinv * expand2To3(p));
	q = shrink3To2(Kinv * expand2To3(q));
	
    A[0] = R0[0][0] - p[0] * R0[0][2];  
    A[1] = R0[1][0] - p[0] * R0[1][2];  
    A[2] = R0[2][0] - p[0] * R0[2][2];
    
    A[3] = R0[0][1] - p[1] * R0[0][2];  
    A[4] = R0[1][1] - p[1] * R0[1][2];  
    A[5] = R0[2][1] - p[1] * R0[2][2];
	
    A[6] = R1[0][0] - q[0] * R1[0][2];  
    A[7] = R1[1][0] - q[0] * R1[1][2];  
    A[8] = R1[2][0] - q[0] * R1[2][2];
	
    A[9] = R1[0][1] - q[1] * R1[0][2];  
    A[10] = R1[1][1] - q[1] * R1[1][2];  
    A[11] = R1[2][1] - q[1] * R1[2][2];
	
    b[0] = t0[2] * p[0] - t0[0];
    b[1] = t0[2] * p[1] - t0[1];
    b[2] = t1[2] * q[0] - t1[0];
    b[3] = t1[2] * q[1] - t1[1];
	
    // Find the least squares result
	int noRowA = 4;
	int noColA = 3;
	int noColB = 1;
	int success = MathUtils::solveLeastSquares(A, b, &noRowA, &noColA, &noColB);
	
	if(success != 0)
	{
		printf("Error! Could not triangulate point!\n");
	}
	
	Vec3f x;
	x[0] = b[0]; x[1] = b[1]; x[2] = b[2];
    
	Vec3f proj0 = R0*x+t0;
	Vec3f proj1 = R1*x+t1;
	
	float dx1 = proj0[0]/proj0[2] - p[0];
	float dy1 = proj0[1]/proj0[2] - p[1];
	float dx2 = proj1[0]/proj1[2] - q[0];
	float dy2 = proj1[1]/proj1[2] - q[1];
	
	error = dx1 * dx1 + dy1 * dy1 + dx2 * dx2 + dy2 * dy2;
	
	return x;
}

bool EpipolarGeometry::isIdentityRotation(Matrix3f R)
{
	double cos_theta = 0.5 * (R[0][0] + R[1][1] + R[2][2] - 1.0);
	cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
	double angle = acos(cos_theta);
	double angleDeg = angle*180.0/M_PI;
	printf("angle:%f\n", angleDeg);
	if(angleDeg < 5.0)
		return true;
	else
		return false;
}

void EpipolarGeometry::findSamplesOnEpipolarLine(Vec3f pixelPos, int pixelIndex, PerspectiveCamera refCam, PerspectiveCamera neighborCam, Img& refImg, Img& neighborImg, int level,
												 float minDepth, float maxDepth, int depthRes, vector<float> &nccEnergies,int minLabelIndex, int totalLabelNo)
{
	int w = neighborImg.width();
	int h = neighborImg.height();
	
	int windowSize = 7;
	//int windowSize = 20 / pow(2.0,level);
	//if(windowSize%2 == 0)
	//	windowSize += 1;
	
	int border = (windowSize-1)/2;
	
	if(pixelPos[0]-border < 0 || pixelPos[0]+border >= refImg.width() || pixelPos[1]-border < 0 || pixelPos[1]+border >= refImg.height())
		return;
	
	Vec3f camCenter = shrink4To3(refCam.getCenter());
	Vec3f rayDir = (refCam.unprojectPixelwithDepth(pixelPos, level) - camCenter);
	
	float delta = (maxDepth - minDepth) / depthRes;
	
	Img refPatch;
	refPatch.resize(windowSize, windowSize);
	for(int x=-border; x<=border; x++)
	{
		for(int y=-border; y<=border; y++)
		{
			refPatch.setColor(x+border, y+border, refImg(pixelPos[0]+x, pixelPos[1]+y));
		}
	}
	
#pragma omp parallel for
	for(int i=0; i<depthRes; i++)
	{
		Vec3f x3d = camCenter + rayDir*(minDepth+i*delta);
		Vec3f x2d = neighborCam.project(expand3To4(x3d), level);
		if(x2d[0] >= border && x2d[0] < w - border && x2d[1] >= border && x2d[1] < h - border)
		{
			Img neighborPatch;
			neighborPatch.resize(windowSize, windowSize);
			for(int x=-border; x<=border; x++)
			{
				for(int y=-border; y<=border; y++)
				{
					neighborPatch.setColor(x+border, y+border, neighborImg(x2d[0]+x, x2d[1]+y));
				}
			}
			
			float score = refPatch.computeNCCScoreWithGrayScale(neighborPatch);
			nccEnergies[pixelIndex*totalLabelNo+minLabelIndex+i] += score;
		}
		else 
		{
			nccEnergies[pixelIndex*totalLabelNo+minLabelIndex+i] += -1.0;
		}
	}
}

Vec3f EpipolarGeometry::computeClosest3DPoint(Vec2f refPixel, Vec2f neighborPixel, PerspectiveCamera &refCam, PerspectiveCamera &neighborCam, int level)
{
	//compute the 3D point that fits best the two image projections
	//minimize: ||P*X - px||^2 where P is the projection matrix, X is the 3D point, px is the image projection
	
	double *A = new double[4*3];
	double *b = new double[4];
	
	//first column
	A[0] = refCam.getElementOfProjectionMatrix(level, 0, 0) - refCam.getElementOfProjectionMatrix(level, 2, 0)*refPixel[0];
	A[1] = refCam.getElementOfProjectionMatrix(level, 1, 0) - refCam.getElementOfProjectionMatrix(level, 2, 0)*refPixel[1];
	A[2] = neighborCam.getElementOfProjectionMatrix(level, 0, 0) - neighborCam.getElementOfProjectionMatrix(level, 2, 0)*neighborPixel[0];
	A[3] = neighborCam.getElementOfProjectionMatrix(level, 1, 0) - neighborCam.getElementOfProjectionMatrix(level, 2, 0)*neighborPixel[1];
	
	//second column
	A[4] = refCam.getElementOfProjectionMatrix(level, 0, 1) - refCam.getElementOfProjectionMatrix(level, 2, 1)*refPixel[0];
	A[5] = refCam.getElementOfProjectionMatrix(level, 1, 1) - refCam.getElementOfProjectionMatrix(level, 2, 1)*refPixel[1];
	A[6] = neighborCam.getElementOfProjectionMatrix(level, 0, 1) - neighborCam.getElementOfProjectionMatrix(level, 2, 1)*neighborPixel[0];
	A[7] = neighborCam.getElementOfProjectionMatrix(level, 1, 1) - neighborCam.getElementOfProjectionMatrix(level, 2, 1)*neighborPixel[1];
	
	//second column
	A[8] = refCam.getElementOfProjectionMatrix(level, 0, 2) - refCam.getElementOfProjectionMatrix(level, 2, 2)*refPixel[0];
	A[9] = refCam.getElementOfProjectionMatrix(level, 1, 2) - refCam.getElementOfProjectionMatrix(level, 2, 2)*refPixel[1];
	A[10] = neighborCam.getElementOfProjectionMatrix(level, 0, 2) - neighborCam.getElementOfProjectionMatrix(level, 2, 2)*neighborPixel[0];
	A[11] = neighborCam.getElementOfProjectionMatrix(level, 1, 2) - neighborCam.getElementOfProjectionMatrix(level, 2, 2)*neighborPixel[1];
	
	
	//b
	b[0] = refCam.getElementOfProjectionMatrix(level, 2, 3)*refPixel[0] - refCam.getElementOfProjectionMatrix(level, 0, 3);
	b[1] = refCam.getElementOfProjectionMatrix(level, 2, 3)*refPixel[1] - refCam.getElementOfProjectionMatrix(level, 1, 3);
	b[2] = neighborCam.getElementOfProjectionMatrix(level, 2, 3)*neighborPixel[0] - neighborCam.getElementOfProjectionMatrix(level, 0, 3);
	b[3] = neighborCam.getElementOfProjectionMatrix(level, 2, 3)*neighborPixel[1] - neighborCam.getElementOfProjectionMatrix(level, 1, 3);
	
	int noRowA = 4;
	int noColA = 3;
	int noColB = 1;
	
	Vec3f point3D(0.0, 0.0, 0.0);
	
	if(MathUtils::solveLeastSquares(A, b, &noRowA, &noColA, &noColB) == 0)
	{
		point3D[0] = b[0];
		point3D[1] = b[1];
		point3D[2] = b[2];
	}
	else 
	{
		printf("EpipolarGeometry: Cannot solve matrix while triangulating point.\n");
	}
	
	delete [] A;
	delete [] b;
	
	return point3D;
}

