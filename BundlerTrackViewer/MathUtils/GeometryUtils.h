/*
 *  GeometryUtils.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _GEOMETRY_UTILS_H
#define _GEOMETRY_UTILS_H

#include "MathUtils.h"
#include "PCA.h"
#include "../BundlerManager/Line.h"

class GeometryUtils
{
public:
	
	static Matrix3f getRotationMatrix(double angle, int axis)
	{
		Matrix3f rot;
		
		switch(axis)
		{
			case 1:
				//x axis
				rot[0][0] = 1;
				rot[1][0] = 0;
				rot[2][0] = 0;
				rot[0][1] = 0; 
				rot[1][1] = cos(angle);
				rot[2][1] = -sin(angle); 
				rot[0][2] = 0;
				rot[1][2] = sin(angle);
				rot[2][2] = cos(angle);
				break;
				
			case 2:
				//y axis
				rot[0][0] = cos(angle);
				rot[1][0] = 0;
				rot[2][0] = sin(angle);
				rot[0][1] = 0; 
				rot[1][1] = 1;
				rot[2][1] = 0; 
				rot[0][2] = -sin(angle);
				rot[1][2] = 0;
				rot[2][2] = cos(angle);
				break;
				
			case 3:
				//z axis
				rot[0][0] = cos(angle);
				rot[1][0] = -sin(angle);
				rot[2][0] = 0;
				rot[0][1] = sin(angle); 
				rot[1][1] = cos(angle);
				rot[2][1] = 0; 
				rot[0][2] = 0;
				rot[1][2] = 0;
				rot[2][2] = 1;
				break;
				
			default:
				rot.setIdentity();
		}
		return rot;
	}
	
	static void computeBoundariesOfTheLine(Vec3f line, int w, int h, Vec2f &minCoord, Vec2f &maxCoord)
	{
		bool intEdge1 = false;
		bool intEdge2 = false;
		bool intEdge3 = false;
		bool intEdge4 = false;
		//edge 1
		float fx1 = 0.0;
		float fy1 = (- line[0] * fx1 - line[2]) / line[1];
		if(0<=fy1 && fy1<=h)
			intEdge1 = true;
		//edge 2
		float fy2 = 0.0;
		float fx2 = (- line[1] * fy2 - line[2]) / line[0];
		if(0<=fx2 && fx2<=w)
			intEdge2 = true;
		//edge 3
		float fx3 = w;
		float fy3 = (- line[0] * fx3 - line[2]) / line[1];
		if(0<=fy3 && fy3<=h)
			intEdge3 = true;
		//edge4
		float fy4 = h;
		float fx4 = (- line[1] * fy4 - line[2]) / line[0];
		if(0<=fx4 && fx4<=w)
			intEdge4 = true;
		
		if(intEdge1 && intEdge2)
		{
			minCoord[0] = fx1;
			minCoord[1] = fy1;
			maxCoord[0] = fx2;
			maxCoord[1] = fy2;
		}
		else if(intEdge1 && intEdge3)
		{
			minCoord[0] = fx1;
			minCoord[1] = fy1;
			maxCoord[0] = fx3;
			maxCoord[1] = fy3;
			
		}
		else if(intEdge1 && intEdge4)
		{
			minCoord[0] = fx1;
			minCoord[1] = fy1;
			maxCoord[0] = fx4;
			maxCoord[1] = fy4;
		}
		else if(intEdge2 && intEdge3)
		{
			minCoord[0] = fx2;
			minCoord[1] = fy2;
			maxCoord[0] = fx3;
			maxCoord[1] = fy3;
		}
		else if(intEdge2 && intEdge4)
		{
			minCoord[0] = fx2;
			minCoord[1] = fy2;
			maxCoord[0] = fx4;
			maxCoord[1] = fy4;
			
		}
		else if(intEdge3 && intEdge4)
		{
			minCoord[0] = fx3;
			minCoord[1] = fy3;
			maxCoord[0] = fx4;
			maxCoord[1] = fy4;
		}	
	}
	
	static void solveForBestFitting2DLine(vector<Line> &lines, vector<float> &weights, Vec3f &bestLine)
	{
		//find an initial guess
		float epsilon = 0.0001;
		int noIter = 5;
		int noLines = lines.size();
		float avgWeight, thresh;
		float totalWeight;
		Vec2f point, avgPoint;
		float length;
		Matrix2f C;
		
		vector<Vec3f> lineEq;
		
		for(int i=0; i<noLines; i++)
		{
			//find normal and offset
			Vec2f lineNormal = lines[i].computeLineNormal2D();
			float lineOffset = lines[i].computeLineDistance2D();
			printf("plot([%f %f], [%f %f], 'r')\n", lines[i].getStartPoint2D()[0], lines[i].getEndPoint2D()[0], lines[i].getStartPoint2D()[1], lines[i].getEndPoint2D()[1]);
			lineEq.push_back(Vec3f(lineNormal[0], lineNormal[1], lineOffset));
		}
		
		for(int i=0; i<noIter; i++)
		{
			totalWeight = 0.0;
			avgPoint[0] = 0.0; avgPoint[1] = 0.0;
			for(int j=0; j<noLines; j++)
			{
				totalWeight += 2.0 * weights[j]*lines[j].computeLineLength2D();
				avgPoint += lines[j].getStartPoint2D()*weights[j]*lines[j].computeLineLength2D();
				avgPoint += lines[j].getEndPoint2D()*weights[j]*lines[j].computeLineLength2D();

			}
			
			if(totalWeight == 0.0)
				return;
			
			avgPoint /= totalWeight;
			
			C.setZero();
			
			//find matrix C
			for(int j=0; j<noLines; j++)
			{
				length = lines[j].computeLineLength2D();
				point = lines[j].getStartPoint2D() - avgPoint;
				C[0][0] = C[0][0] + length*point[0]*point[0];
				C[0][1] = C[0][1] + length*point[0]*point[1];
				C[1][0] = C[1][0] + length*point[1]*point[0];
				C[1][1] = C[1][1] + length*point[1]*point[1];
					
				point = lines[j].getEndPoint2D() - avgPoint;
				C[0][0] = C[0][0] + length*point[0]*point[0];
				C[0][1] = C[0][1] + length*point[0]*point[1];
				C[1][0] = C[1][0] + length*point[1]*point[0];
				C[1][1] = C[1][1] + length*point[1]*point[1];
			}
			
			Vec2f eigenValues;
			Matrix2f eigenVectors;
			int noRows = 2;
			
			MathUtils::eigenValueAnalysis(&(C[0][0]), &noRows, &(eigenValues[0]), &(eigenVectors[0][0]));
			bestLine[0] = eigenVectors[1][0];
			bestLine[1] = eigenVectors[1][1];
			bestLine[2] = 0.0 - bestLine[0]*avgPoint[0] - bestLine[1]*avgPoint[1];
			
			Vec2f dir(bestLine[1], -bestLine[0]);
			
			avgWeight = 0.0;
			for(int j=0; j<noLines; j++)
			{
				float length1 = abs(dot(Vec2f(bestLine[0], bestLine[1]), lines[j].getStartPoint2D()) + bestLine[2]);
				float length2 = abs(dot(Vec2f(bestLine[0], bestLine[1]), lines[j].getEndPoint2D()) + bestLine[2]);
				float dist = (length1+length2)/2.0;
				//float area = findAreaUnderLine(Vec2f(bestLine[0], bestLine[1]), bestLine[2], lines[j].getStartPoint2D(), lines[j].getEndPoint2D());
				weights[j] = 1.0 / (dist+epsilon);
				//weights[j] = 1.0 / (area+epsilon);
				//weights[j] = dist;
				avgWeight += weights[j];
				//printf("dist:%f\n", dist);
			}
			avgWeight /= noLines;
			thresh = avgWeight / 10.0;
			for(int j=0; j<noLines; j++)
			{
				if(weights[j] < thresh)
					weights[j] = 0.0;
			}
		}
	}
	
	static void findIntersectionOfContourLines(vector<Vec3f> &contourLines, vector<Vec2f> &intersectionPts)
	{
		int noLines = contourLines.size();
		
		for(int i=0; i<noLines; i++)
		{
			int j=i+1;
			if(i==noLines-1)
				j=0;
			
			printf("line1:%f,%f,%f line2:%f,%f,%f\n", contourLines[i][0], contourLines[i][1], contourLines[i][2],
													  contourLines[j][0], contourLines[j][1], contourLines[j][2]);
			//find intersection between line i and i+1
			float denom = contourLines[i][0]*contourLines[j][1] - contourLines[j][0]*contourLines[i][1];
			if(denom !=  0.0)
			{
				float x = (contourLines[j][1]*(-contourLines[i][2])-contourLines[i][1]*(-contourLines[j][2])) / denom;
				float y = (contourLines[i][0]*(-contourLines[j][2])-contourLines[j][0]*(-contourLines[i][2])) / denom;
				intersectionPts.push_back(Vec2f(x,y));
				printf("intersection:%f %f\n", x,y);
			}
			
		}
	}
	
	static bool clipLineWrt2DRectangle(Vec3f line, float xmin, float xmax, float ymin, float ymax, Vec2f &minCoord, Vec2f &maxCoord)
	{
		bool intEdge1 = false;
		bool intEdge2 = false;
		bool intEdge3 = false;
		bool intEdge4 = false;
		minCoord[0] = 0.0; minCoord[1] = 0.0;
		maxCoord[0] = 0.0; maxCoord[1] = 0.0;
		//edge 1
		float fx1 = xmin;
		float fy1 = (- line[0] * fx1 - line[2]) / line[1];
		if(ymin<=fy1 && fy1<=ymax)
			intEdge1 = true;
		//edge 2
		float fy2 = ymin;
		float fx2 = (- line[1] * fy2 - line[2]) / line[0];
		if(xmin<=fx2 && fx2<=xmax)
			intEdge2 = true;
		//edge 3
		float fx3 = xmax;
		float fy3 = (- line[0] * fx3 - line[2]) / line[1];
		if(ymin<=fy3 && fy3<=ymax)
			intEdge3 = true;
		//edge4
		float fy4 = ymax;
		float fx4 = (- line[1] * fy4 - line[2]) / line[0];
		if(xmin<=fx4 && fx4<=xmax)
			intEdge4 = true;
		
		if(intEdge1 && intEdge2)
		{
			minCoord[0] = fx1;
			minCoord[1] = fy1;
			maxCoord[0] = fx2;
			maxCoord[1] = fy2;
			return true;
		}
		else if(intEdge1 && intEdge3)
		{
			minCoord[0] = fx1;
			minCoord[1] = fy1;
			maxCoord[0] = fx3;
			maxCoord[1] = fy3;
			return true;
		}
		else if(intEdge1 && intEdge4)
		{
			minCoord[0] = fx1;
			minCoord[1] = fy1;
			maxCoord[0] = fx4;
			maxCoord[1] = fy4;
			return true;
		}
		else if(intEdge2 && intEdge3)
		{
			minCoord[0] = fx2;
			minCoord[1] = fy2;
			maxCoord[0] = fx3;
			maxCoord[1] = fy3;
			return true;
		}
		else if(intEdge2 && intEdge4)
		{
			minCoord[0] = fx2;
			minCoord[1] = fy2;
			maxCoord[0] = fx4;
			maxCoord[1] = fy4;
			return true;
		}
		else if(intEdge3 && intEdge4)
		{
			minCoord[0] = fx3;
			minCoord[1] = fy3;
			maxCoord[0] = fx4;
			maxCoord[1] = fy4;
			return true;
		}
		else
		{
			return false;
		}
	}
	
	static bool isPointInsidePolygonIn2D(vector<Vec2f> &polygonCoord, Vec2f point)
	{
		int noVertices = polygonCoord.size();
		bool isInside = false;
		int i;
		int j = noVertices - 1;
		for(i=0; i<noVertices; i++)
		{
			if(floor(polygonCoord[i][0]+0.5) == floor(point[0]+0.5) && floor(polygonCoord[i][1]+0.5) == floor(point[1]+0.5))
				return true;
			
			if ((polygonCoord[i][1]<point[1] && polygonCoord[j][1]>=point[1]) || (polygonCoord[j][1]<point[1] && polygonCoord[i][1]>=point[1]))
			{
				if(polygonCoord[i][0]<=point[0] || polygonCoord[j][0]<=point[0])
				{
					//intersection point
					float intX = polygonCoord[i][0] + (point[1]-polygonCoord[i][1])*(polygonCoord[j][0]-polygonCoord[i][0])/(polygonCoord[j][1]-polygonCoord[i][1]); 
					if (intX < point[0]) 
					{
						isInside=!isInside;
					}
				}
			}
			j=i;
		}
		
		if(!isInside)
		{
			//printf("p:%f %f\n", point[0], point[1]);
		}
		return isInside;
		
	}
	
	static void findRigidTransformation(vector<Vec3f> &sourcePoints, vector<Vec3f> &targetPoints, Matrix3f &rot, Vec3f &trans)
	{
		int noSource = sourcePoints.size();
		int noTarget = targetPoints.size();
		
		Vec3f avgSource(0.0, 0.0, 0.0);
		Vec3f avgTarget(0.0, 0.0, 0.0);
		
		for(int i=0; i<noSource; i++)
		{
			avgSource += sourcePoints[i];
		}
		avgSource /= noSource;
		
		for(int i=0; i<noTarget; i++)
		{
			avgTarget += targetPoints[i];
		}
		avgTarget /= noTarget;
		
		float *XMat = new float[noSource*3];
		float *YMatTrans = new float[noTarget*3];
		
		for(int i=0; i<noSource; i++)
		{
			Vec3f p = sourcePoints[i];
			p = p - avgSource;
			XMat[i] = p[0];
			XMat[noSource + i] = p[1];
			XMat[noSource*2 + i] = p[2];
		}
		
		//transpose
		for(int i=0; i<noTarget; i++)
		{
			Vec3f p = targetPoints[i];
			p = p - avgTarget;
			YMatTrans[i*3] = p[0];
			YMatTrans[i*3 + 1] = p[1];
			YMatTrans[i*3 + 2] = p[2];
		}
		
		//X*Y_trans
		double *mat = new double[9];
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				double val = 0.0;
				for(int k=0; k<noSource; k++)
				{
					val += XMat[i*noSource+k]*YMatTrans[k*3+j];
				}
				mat[j*3+i] = val;
			}
		}
		
		int noRows = 3;
		int noCols = 3;
		double *u = new double[9];
		double *s = new double[9];
		double *vt = new double[9];
		MathUtils::computeSVD(mat, &noRows, &noCols, u, s, vt);
		
		//R = v*uT
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				float val = 0.0;
				for(int k=0; k<3; k++)
				{
					val += vt[i*3+k]*u[k*3+j];
				}
				rot[j][i] = val;
			}
		}
		
		printf("det:%f\n", rot.getDeterminant());
		
		/*if(rot.getDeterminant() == -1.0)
		{
			u[6] = -u[6];
			u[7] = -u[7];
			u[8] = -u[8];
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; j++)
				{
					float val = 0.0;
					for(int k=0; k<3; k++)
					{
						val += vt[i*3+k]*u[k*3+j];
					}
					rot[j][i] = val;
				}
			}
		}*/
		
		trans = avgTarget - rot*avgSource;
	}
	
	static void projectDataToPlane(vector<Vec3f> &points, Vec3f &n, Vec3f &centroid)
	{
		//fit a plane to the data
		int noPoints = points.size();
		
		PCA3f pca;
		Vec3f avgNormal(0.0,0.0,0.0);
		centroid = Vec3f(0.0, 0.0, 0.0);
		for(int i=0; i<noPoints; i++)
		{
			pca.addPoint(points[i]);
		}
		
		//analyze eigen structure
		Vec3f eigenValues;
		Matrix3f eigenVectors;
		pca.analyze(eigenValues, eigenVectors, centroid);
		
		n = eigenVectors[2];
		n = n.normalize();
		
		for(int i=0; i<noPoints; i++)
		{
			Vec3f pos = points[i];
			points[i] = centroid + ((pos-centroid) - n*dot(n, pos-centroid));
		}
	}
	
};

#endif