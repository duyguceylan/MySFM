/*
 *  Plane.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 6/29/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _PLANE_H
#define _PLANE_H

#include "../MathUtils/MyMatrix.h"

class Plane
{
protected:
	Vec3f normal;
	float distance;
	Vec3f xAxis, yAxis;
	Matrix3f rotMatrix;
	Vec3f center;
	
	Vec3f centroid;
	Vec3f minCoord;
	Vec3f maxCoord;
	
public:
	Plane() {}
	Plane(Vec3f normal, float distance);
	
	void flipNormal() {normal = -normal; distance = -distance; }
	void setAxis(Vec3f xAxis_, Vec3f yAxis_);
	void computeAxis();
	void setPlaneCentroid(Vec3f centroid_) {centroid = centroid_;}
	
	Vec3f getPlaneNormal() {return normal; }
	float getPlaneDistance() {return distance; }
	Vec3f getPlaneXAxis() {return xAxis; }
	Vec3f getPlaneYAxis() {return yAxis; }
	Vec3f getPlaneCentroid() {return centroid; }
	
	bool findPlaneRayIntersection(Vec3f rayOrg, Vec3f rayDir, Vec3f &intersection);
	bool findPlaneRayIntersectionWithTranslation(Vec3f rayOrg, Vec3f rayDir, float translation, Vec3f &intersection);
	Vec3f findProjectionOnPlane(Vec3f point);
	Vec2f transformPointToLocalCoordinateSystem(Vec3f pt3D);
	Vec3f transformPointToGlobalCoordinateSystem(Vec2f pt2D);

	void setBoundingBox(Vec3f minCoord_, Vec3f maxCoord_) {minCoord = minCoord_; maxCoord = maxCoord_;}
	void getBoundingBox(Vec3f &minCoord_, Vec3f &maxCoord_) {minCoord_ = minCoord; maxCoord_ = maxCoord;}
};

#endif