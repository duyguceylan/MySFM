/*
 *  PointCloud.h
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/1/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _POINT_CLOUD_
#define _POINT_CLOUD_

#include "Data3D.h"

#include <iostream>
#include <vector>

typedef OpenMesh::Vec2i Vec2i;

class PointCloud : public Data3D
{
private:	
	
public:
	PointCloud();
	PointCloud(const char *filename);
	~PointCloud();
	
	void updateVertexNormals(int noNeighbors=20, bool usePreviousOrientation=true);
	
	PointCloud* getTranslatedData(Vec3f trans);
	PointCloud* getTransformedData(Matrix4f mat);
	
	PointCloud& operator=(const PointCloud &rhs);
};

#endif