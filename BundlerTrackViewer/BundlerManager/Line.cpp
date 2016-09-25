/*
 *  Line.cpp
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 6/29/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include "Line.h"
#include "../MathUtils/MathUtils.h"

Line::Line()
{
	matched3DLine = false;
	isContour = false;
	matches.clear();
	startPoint3D = Vec3f(0.0, 0.0, 0.0);
	endPoint3D = Vec3f(0.0, 0.0, 0.0);
	startPoint2D = Vec2f(0.0, 0.0);
	endPoint2D = Vec2f(0.0, 0.0);
	clippedStartPoint = Vec2f(0.0, 0.0);
	clippedEndPoint = Vec2f(0.0, 0.0);
}

Line::Line(Vec3f startPoint_, Vec3f endPoint_)
{
	startPoint3D = startPoint_;
	endPoint3D = endPoint_;
	matched3DLine = false;
	isContour = false;
	matches.clear();
}

Line::Line(Vec2f projStartPoint_, Vec2f projEndPoint_)
{
	startPoint2D = projStartPoint_;
	endPoint2D = projEndPoint_;
	clippedStartPoint = startPoint2D;
	clippedEndPoint = endPoint2D;
	matched3DLine = false;
	isContour = false;
	matches.clear();
}

void Line::projectToImage(Matrix4f projMatrix)
{
	Vec3f p = shrink4To3(projMatrix * expand3To4(startPoint3D));
	startPoint2D = Vec2f(0.0, 0.0);
	if(p[2] != 0.0)
	{
		startPoint2D[0] = p[0]/p[2];
		startPoint2D[1] = p[1]/p[2];
	}

	p = shrink4To3(projMatrix * expand3To4(endPoint3D));
	endPoint2D = Vec2f(0.0, 0.0);
	if(p[2] != 0.0)
	{
		endPoint2D[0] = p[0]/p[2];
		endPoint2D[1] = p[1]/p[2];
	}
}

int Line::doesImageHaveProjection(int imgIndex)
{
	for(int i=0; i<imageProjections.size(); i++)
	{
		if(imageProjections[i].imgIndex == imgIndex)
		{
			return i;
		}
	}
	return -1;
}

Vec2f Line::generateRandomPointOnLine2D()
{
	Vec2f dir = endPoint2D - startPoint2D;
	return startPoint2D + dir*((float)rand()/(float)RAND_MAX);
}

float Line::computeLineLength2D()
{
	return (endPoint2D - startPoint2D).length();
}

Vec2f Line::computeLineNormal2D()
{
	Vec2f n;
	n[0] = startPoint2D[1] - endPoint2D[1];
	n[1] = endPoint2D[0] - startPoint2D[0];
	n.normalize();
	return n;
}

float Line::computeLineDistance2D()
{
	Vec2f n = computeLineNormal2D();
	return 0.0 - n[0]*startPoint2D[0] - n[1]*startPoint2D[1];
}

float Line::findPerpendicularDistanceTo2DLine(Vec2f pixelPos)
{
	Vec2f n = computeLineNormal2D();
	Vec2f r = pixelPos - startPoint2D;
	return dot(n, r);
}

void Line::findBoundingRectangle(vector<Vec2f> &boundingRect)
{
	float threshold = 4.0;
    
    Vec2f lineDir = (endPoint2D - startPoint2D).normalize();
	Vec2f lineNormal = computeLineNormal2D();
    
    boundingRect.push_back(startPoint2D - lineDir*threshold + lineNormal*threshold);
	boundingRect.push_back(endPoint2D + lineDir*threshold + lineNormal*threshold);
	boundingRect.push_back(endPoint2D + lineDir*threshold - lineNormal*threshold);
	boundingRect.push_back(startPoint2D - lineDir*threshold - lineNormal*threshold);
}

bool Line::isInsideBoundingRectangle(vector<Vec2f> &boundingRect, Vec2f pt)
{
	int fp, fp1, fp2, fp3, fp4; 
	//edge 1
	if((boundingRect[1][0]-boundingRect[0][0]) == 0)
	{
        fp = pt[0] - boundingRect[0][0];
		fp3 = boundingRect[2][0] - boundingRect[0][0];
	}
    else
	{
        fp = pt[1] - boundingRect[0][1] - ((boundingRect[1][1]-boundingRect[0][1]) / (boundingRect[1][0]-boundingRect[0][0]))*(pt[0] - boundingRect[0][0]);
		fp3 = boundingRect[2][1] - boundingRect[0][1] - ( (boundingRect[1][1]-boundingRect[0][1]) / (boundingRect[1][0]-boundingRect[0][0]) )*(boundingRect[2][0] - boundingRect[0][0]);
	}
    
	if(fp3*fp < 0)
		return false;
    
    //edge 2
	if((boundingRect[2][0]-boundingRect[1][0]) == 0)
	{
        fp = pt[0] - boundingRect[1][0];
		fp4 = boundingRect[3][0] - boundingRect[1][0];
	}
    else
	{
        fp = pt[1] - boundingRect[1][1] - ( (boundingRect[2][1]-boundingRect[1][1]) / (boundingRect[2][0]-boundingRect[1][0]) )*(pt[0] - boundingRect[1][0]);
		fp4 = boundingRect[3][1] - boundingRect[1][1] - ( (boundingRect[2][1]-boundingRect[1][1]) / (boundingRect[2][0]-boundingRect[1][0]) )*(boundingRect[3][0] - boundingRect[1][0]);
	}
    if(fp4*fp < 0)
		return false;
   
	//edge 3
    if((boundingRect[3][0]-boundingRect[2][0]) == 0)
	{
        fp = pt[0] - boundingRect[2][0];
		fp1 = boundingRect[0][0] - boundingRect[2][0];
	}
    else
	{
        fp = pt[1] - boundingRect[2][1] - ( (boundingRect[3][1]-boundingRect[2][1]) / (boundingRect[3][0]-boundingRect[2][0]) )*(pt[0] - boundingRect[2][0]);
		fp1 = boundingRect[0][1] - boundingRect[2][1] - ( (boundingRect[3][1]-boundingRect[2][1]) / (boundingRect[3][0]-boundingRect[2][0]) )*(boundingRect[0][0] - boundingRect[2][0]);
	}
    if(fp1*fp < 0)
		return false;
   
	//edge 4
    if((boundingRect[0][0]-boundingRect[3][0]) == 0)
	{
        fp = pt[0] - boundingRect[3][0];
		fp2 = boundingRect[1][0] - boundingRect[3][0];
	}
    else
	{
        fp = pt[1] - boundingRect[3][1] - ( (boundingRect[0][1]-boundingRect[3][1]) / (boundingRect[0][0]-boundingRect[3][0]) )*(pt[0] - boundingRect[3][0]);
		fp2 = boundingRect[1][1] - boundingRect[3][1] - ( (boundingRect[0][1]-boundingRect[3][1]) / (boundingRect[0][0]-boundingRect[3][0]) )*(boundingRect[1][0] - boundingRect[3][0]);
	}
    if(fp2*fp < 0)
		return false;
    
	return true;
}

Vec3f Line::compute2DLineNormalFormParams()
{
	Vec3f p;
	if(startPoint2D[0] == endPoint2D[0])
	{
		if(startPoint2D[0]>=0.0)
			p[0] = 1.0;
		else 
			p[0] = -1.0;
		p[1] = 0.0;
		p[2] = 0.0 - p[0] * startPoint2D[0];
	}
	else if(startPoint2D[1] == endPoint2D[1])
	{
		p[0] = 0.0;
		if(startPoint2D[1]>=0.0)
			p[1] = 1.0;
		else
			p[1] = -1.0;
		p[2] = 0.0 - p[1]*startPoint2D[1];
	}
	else 
	{
		//slope - intercept
		float m = (endPoint2D[1] - startPoint2D[1]) / (endPoint2D[0] - startPoint2D[0]);
		float b = startPoint2D[1] - m*startPoint2D[0];
		Vec2f yInt(0.0, b);
		Vec2f xInt(-b/m, 0.0);
		float length = sqrt(yInt[1]*yInt[1] + xInt[0]*xInt[0]);
		float cosAlpha = yInt[1] / length;
		float sinAlpha = xInt[0] / length;
		if((yInt[1]>0.0 && xInt[0]<0.0) || (yInt[1]<0.0 && xInt[0]>0.0))
		{
			cosAlpha *= -1.0;
			sinAlpha *= -1.0;
		}
		p[0] = cosAlpha;
		p[1] = sinAlpha;
		p[2] = 0.0 - p[0]*startPoint2D[0] - p[1]*startPoint2D[1];
	}
	return p;
}

//return
//0:parallel
//1: intersecting
//2:non-coplanar
LINERELATION Line::findRelationWith3DLine(Line *line, Vec3f &intersectionPt)
{
	float epsilon = 0.5;
	float angleThreshold = 20.0*M_PI/180.0;
	
	//find the line directions
	Vec3f lineDir1= endPoint3D - startPoint3D;
	Vec3f lineDir2 = line->getEndPoint3D() - line->getStartPoint3D();
	
	//line0: bo+s*mo   line1: b1+s*m1
	float a = dot(lineDir1, lineDir1);
	float b = -dot(lineDir1, lineDir2);
	float c = dot(lineDir2, lineDir2);
	float d = dot(lineDir1, (startPoint3D - line->getStartPoint3D()));
	float e = -dot(lineDir2, (startPoint3D - line->getStartPoint3D()));
	float f = dot((startPoint3D - line->getStartPoint3D()), (startPoint3D - line->getStartPoint3D()));
	
	float det = a*c - b*b;
	Vec3f dir1 = lineDir1;
	Vec3f dir2 = lineDir2;
	if(abs(dot(dir1.normalize(), dir2.normalize())) > cos(angleThreshold))
	{
		//line are same?
		Vec3f projPoint = startPoint3D + dir1 * dot(dir1, line->getStartPoint3D() - startPoint3D);
		if((projPoint-line->getStartPoint3D()).length() < epsilon)
		   return SAME;
		
		//lines are parallel
		else
		   return PARALLEL;
		
	}
	
	//s = b*e - c*d / det   t = b*d - a*e / det
	float s = b*e - c*d / det;
	float t = b*d - a*e / det;
	
	float dist = (startPoint3D + lineDir1*s - line->getStartPoint3D() - lineDir2*t).length();
	//printf("dist:%f\n", dist);
	if(dist < epsilon)
	{
		intersectionPt = startPoint3D + s*lineDir1;
		return INTERSECTING;
	}
	else 
	{
		return NONCOPLANAR;
	}
}

LINERELATION Line::findRelationWith2DLine(Line *line, Vec2f &intersectionPt, bool extend)
{
	Vec2f neighStart = line->getStartPoint2D();
	Vec2f neighEnd = line->getEndPoint2D();
	
	float A1 = (endPoint2D - startPoint2D)[1];
    float B1 = (startPoint2D-endPoint2D)[0];
    float C1 = A1*endPoint2D[0]+B1*endPoint2D[1];
    
    float A2 = (neighEnd-neighStart)[1];
    float B2 = (neighStart-neighEnd)[0];
    float C2 = A2*neighEnd[0]+B2*neighEnd[1];
    
    float det = A1*B2 - A2*B1;
    
	float epsilon = (endPoint2D - startPoint2D).length() / 100;
	
    if(abs(det) < 1E-3)
		return PARALLEL;
    
    intersectionPt[0] = (B2*C1 - B1*C2)/det;
	intersectionPt[1] = (A1*C2 - A2*C1)/det;
	
	printf("intersection:%f %f\n", intersectionPt[0], intersectionPt[1]);
	
	if(!extend)
	{
		if(((startPoint2D[0]<intersectionPt[0] || abs(startPoint2D[0]-intersectionPt[0])<epsilon) && (intersectionPt[0]<endPoint2D[0] || abs(intersectionPt[0]-endPoint2D[0])<epsilon)) 
		   || ((endPoint2D[0]<intersectionPt[0] ||  abs(endPoint2D[0]<intersectionPt[0])<epsilon) && (intersectionPt[0]<startPoint2D[0] || abs(intersectionPt[0]-startPoint2D[0])<epsilon)))
		{
			if(((startPoint2D[1]<intersectionPt[1] || abs(startPoint2D[1]-intersectionPt[1])<epsilon) && (intersectionPt[1]<endPoint2D[1] || abs(intersectionPt[1]-endPoint2D[1])<epsilon)) 
			   || ((endPoint2D[1]<intersectionPt[1] || abs(endPoint2D[1]-intersectionPt[1])<epsilon) && (intersectionPt[1]<startPoint2D[1] || abs(intersectionPt[1]-startPoint2D[1])<epsilon)))
				{
					if(((neighStart[0]<intersectionPt[0] || abs(neighStart[0]-intersectionPt[0])<epsilon) && (intersectionPt[0]<neighEnd[0] || abs(intersectionPt[0]-neighEnd[0])<epsilon)) 
					   || ((neighEnd[0]<intersectionPt[0] || abs(neighEnd[0]-intersectionPt[0])<epsilon) && (intersectionPt[0]-neighStart[0] || abs(intersectionPt[0]-neighStart[0])<epsilon)))
					{
						if(((neighStart[1]<intersectionPt[1] || abs(neighStart[1]-intersectionPt[1])<epsilon) && (intersectionPt[1]<neighEnd[1] || abs(intersectionPt[1]-neighEnd[1])<epsilon)) || 
						   ((neighEnd[1]<intersectionPt[1] || abs(neighEnd[1]-intersectionPt[1])<epsilon) && (intersectionPt[1]<neighStart[1] || abs(intersectionPt[1]-neighStart[1])<epsilon)))
						{
							return INTERSECTING;
						}
						else
							return NONINTERSECTING;
					}
					else
						return NONINTERSECTING;
				}
			else
				return NONINTERSECTING;
		}
		else
			return NONINTERSECTING;
	}
	else
	{
		return INTERSECTING;
	}

}

bool Line::isVisibleByCamera(Matrix4f projMatrix, Vec2f &projStart, Vec2f &projEnd, int width, int height)
{
	bool start, end;
	//check start point
	Vec3f proj = shrink4To3(projMatrix * expand3To4(startPoint3D));
	if(proj[2] == 0.0)
		return false;
	proj[0] /= proj[2];
	proj[1] /= proj[2];
	projStart = shrink3To2(proj);
	
	printf("start:%f %f\n", proj[0], proj[1]);
	
	if(proj[0]<0.0 || proj[0]>=width || proj[1]<0.0 || proj[1]>=height)
		start = false;
	else
		start = true;
	
	//check end point
	proj = shrink4To3(projMatrix * expand3To4(endPoint3D));
	if(proj[2] == 0.0)
		return false;
	proj[0] /= proj[2];
	proj[1] /= proj[2];
	projEnd = shrink3To2(proj);
	
	printf("end:%f %f\n", proj[0], proj[1]);
	
	if(proj[0]<0.0 || proj[0]>=width || proj[1]<0.0 || proj[1]>=height)
		end = false;
	else
		end = true;
	
	if(start && end)
		return true;
	
	bool intEdge1 = false;
	bool intEdge2 = false;
	bool intEdge3 = false;
	bool intEdge4 = false;
	//edge 1
	float fx1 = 0.0;
	float t1 = (- projStart[0])/(projEnd[0]-projStart[0]);
	float fy1 = projStart[1] + t1*(projEnd[1]-projStart[1]);
	if(0<=fy1 && fy1<=height && t1>=0.0 && t1<=1.0)
		intEdge1 = true;
	//edge 2
	float fy2 = 0.0;
	float t2 = (- projStart[1])/(projEnd[1]-projStart[1]);
	float fx2 = projStart[0] + t2*(projEnd[0]-projStart[0]);
	if(0<=fx2 && fx2<=width && t2>=0.0 && t2<=1.0)
		intEdge2 = true;
	//edge 3
	float fx3 = width;
	float t3 = (width - projStart[0])/(projEnd[0]-projStart[0]);
	float fy3 = projStart[1] + t3*(projEnd[1]-projStart[1]);
	if(0<=fy3 && fy3<=height && t3>=0.0 && t3<=1.0)
		intEdge3 = true;
	//edge4
	float fy4 = height;
	float t4 = (height - projStart[1])/(projEnd[1]-projStart[1]);
	float fx4 = projStart[0] + t4*(projEnd[0]-projStart[0]);
	if(0<=fx4 && fx4<=width && t4>=0.0 && t4<=1.0)
		intEdge4 = true;
	
	if(intEdge1 && intEdge2)
	{
		projStart[0] = fx1;
		projStart[1] = fy1;
		projEnd[0] = fx2;
		projEnd[1] = fy2;
		return true;
	}
	else if(intEdge1 && intEdge3)
	{
		projStart[0] = fx1;
		projStart[1] = fy1;
		projEnd[0] = fx3;
		projEnd[1] = fy3;
		return true;
	}
	else if(intEdge1 && intEdge4)
	{
		projStart[0] = fx1;
		projStart[1] = fy1;
		projEnd[0] = fx4;
		projEnd[1] = fy4;
		return true;
	}
	else if(intEdge2 && intEdge3)
	{
		projStart[0] = fx2;
		projStart[1] = fy2;
		projEnd[0] = fx3;
		projEnd[1] = fy3;
		return true;
	}
	else if(intEdge2 && intEdge4)
	{
		projStart[0] = fx2;
		projStart[1] = fy2;
		projEnd[0] = fx4;
		projEnd[1] = fy4;
		return true;
	}
	else if(intEdge3 && intEdge4)
	{
		projStart[0] = fx3;
		projStart[1] = fy3;
		projEnd[0] = fx4;
		projEnd[1] = fy4;
		return true;
	}
	else if(intEdge1)
	{
		if(start && !end)
		{
			projEnd[0] = fx1;
			projEnd[1] = fy1;
			return true;
		}
		else if(end && !start)
		{
			projStart[0] = fx1;
			projStart[1] = fy1;
			return true;
		}
	}
	else if(intEdge2)
	{
		if(start && !end)
		{
			projEnd[0] = fx2;
			projEnd[1] = fy2;
			return true;
		}
		else if(end && !start)
		{
			projStart[0] = fx2;
			projStart[1] = fy2;
			return true;
		}
	}
	else if(intEdge3)
	{
		if(start && !end)
		{
			projEnd[0] = fx3;
			projEnd[1] = fx3;
			return true;
		}
		else if(end && !start)
		{
			projStart[0] = fx3;
			projStart[1] = fy3;
			return true;
		}
	}
	else if(intEdge4)
	{
		if(start && !end)
		{
			projEnd[0] = fx4;
			projEnd[1] = fy4;
			return true;
		}
		else if(end && !start)
		{
			projStart[0] = fx4;
			projStart[1] = fy4;
			return true;
		}
	}
	else 
	{
		return false;
	}
}

void Line::fillLineMask(int width, int height, Img *tmp)
{
	//find slope
	float yDiff = abs(endPoint2D[1] - startPoint2D[1]);
	float xDiff = abs(endPoint2D[0] - startPoint2D[0]);
	float sx, sy;
	if(startPoint2D[0] < endPoint2D[0])
		sx = 1;
	else
		sx = -1;
	if(startPoint2D[1] < endPoint2D[1])
		sy = 1;
	else
		sy = -1;
	
	float error = xDiff - yDiff;
	int x = startPoint2D[0];
	int y = startPoint2D[1];
	
	while(true)
	{
		if(x<0 || x>=tmp->width() || y<0 || y>=tmp->height())
			break;
		
		tmp->setColor(x,y, Color(1.0,0.0,0.0));
		if(x == (int)endPoint2D[0] && y==(int)endPoint2D[1])
			break;
		else if(sx==1 && sy==1 && x >=endPoint2D[0] && y >=endPoint2D[1])
			break;
		else if(sx==-1 && sy==1 && x <=endPoint2D[0] && y >=endPoint2D[1])
			break;
		else if(sx==1 && sy==-1 && x >=endPoint2D[0] && y<=endPoint2D[1])
			break;
		else if(sx==-1 && sy==-1 && x <=endPoint2D[0] && y <=endPoint2D[1])
			break;
		
		float e2 = 2.0*error;
		if (e2 > -yDiff)
		{
			error -= yDiff;
			x += sx;
		}
		if (e2 < xDiff)
		{
			error += xDiff;
			y += sy ;
		}
	}
}