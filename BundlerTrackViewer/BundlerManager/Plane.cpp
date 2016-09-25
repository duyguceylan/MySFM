/*
 *  Plane.cpp
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 6/29/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "Plane.h"

Plane::Plane(Vec3f normal_, float distance_)
{
	normal = normal_;
	distance = distance_;
	center = Vec3f(0.0, 0.0, 0.0);
	rotMatrix.setZero();
}

void Plane::computeAxis()
{
	Vec3f yAxis_(0.0, 1.0, 0.0);
	Vec3f xAxis_ = cross(yAxis_, normal).normalize();
	yAxis_ = cross(normal, xAxis_).normalize();
	xAxis_ = cross(yAxis_, normal).normalize();
	setAxis(xAxis_, yAxis_);
}

void Plane::setAxis(Vec3f xAxis_, Vec3f yAxis_)
{
	xAxis = xAxis_;
	yAxis = yAxis_;
	rotMatrix[0] = xAxis;
	rotMatrix[1] = yAxis;
	rotMatrix[2] = normal;
	rotMatrix = rotMatrix.transpose();
	float epsilon = 1E-5;
	if(fabs(normal[2]) > epsilon)
		center = Vec3f(0.0, 0.0, -distance/normal[2]);
	else if(fabs(normal[1]) > epsilon)
		center = Vec3f(0.0, -distance/normal[1], 0.0);
	else if(fabs(normal[0]) > epsilon)
		center = Vec3f(-distance/normal[0], 0.0, 0.0);

}

bool Plane::findPlaneRayIntersection(Vec3f rayOrg, Vec3f rayDir, Vec3f &intersection)
{
	//intersection point: rayOrg + t*rayDir
	//normal^T(rayOrg + t*rayDir) + distance = 0.0
	
	//check if there's an intersection
	if(dot(normal, rayDir) == 0.0)
		return false;
		
	float t = (0.0 - distance - dot(normal, rayOrg)) / (dot(normal, rayDir));
	intersection = rayOrg + rayDir * t;
	return true;
}

bool Plane::findPlaneRayIntersectionWithTranslation(Vec3f rayOrg, Vec3f rayDir, float translation, Vec3f &intersection)
{
	float newDistance = distance - dot(normal, normal*translation);
	//printf("newDistance:%f\n", newDistance);
	//check if there's an intersection
	if(dot(normal, rayDir) == 0.0)
		return false;
	
	float t = (0.0 - newDistance - dot(normal, rayOrg)) / (dot(normal, rayDir));
	intersection = rayOrg + rayDir * t;
	return true;
}

Vec3f Plane::findProjectionOnPlane(Vec3f point)
{
	float dist = distance + dot(normal, point);
	return point - normal*dist;
}

Vec2f Plane::transformPointToLocalCoordinateSystem(Vec3f pt3D)
{
	return shrink3To2(rotMatrix * (pt3D - center));
}

Vec3f Plane::transformPointToGlobalCoordinateSystem(Vec2f pt2D)
{
	return rotMatrix.transpose()*Vec3f(pt2D[0], pt2D[1], 0.0) + center;
}

