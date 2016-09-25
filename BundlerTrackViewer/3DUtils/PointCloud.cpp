/*
 *  PointCloud.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 3/13/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "PointCloud.h"

PointCloud::PointCloud() : Data3D() 
{
	dataType = POINTCLOUD;
}

PointCloud::PointCloud(const char *filename) : Data3D(filename) 
{
	dataType = POINTCLOUD;
}

PointCloud::~PointCloud()
{
}

PointCloud& PointCloud::operator=(const PointCloud &rhs)
{
	if(this != &rhs)
	{
		meshData = rhs.meshData;
		dataType = rhs.dataType;
		modelName = rhs.modelName;
	}
	
	return *this;	
}

PointCloud* PointCloud::getTransformedData(Matrix4f mat)
{
	PointCloud *result = new PointCloud();
	
	Vec3f pos, normal;
	
	for(int i=0; i<meshData.n_vertices(); i++)
	{
		pos = shrink4To3(mat * expand3To4(meshData.point(meshData.vertex_handle(i))));
		normal = shrink4To3(mat) * meshData.normal(meshData.vertex_handle(i));
		result->addPoint(pos, normal);
	}
	return result;
}

PointCloud* PointCloud::getTranslatedData(Vec3f trans)
{
	PointCloud *result = new PointCloud();
	
	for(int i=0; i<meshData.n_vertices(); i++)
	{
		result->addPoint(meshData.point(meshData.vertex_handle(i))+trans, meshData.normal(meshData.vertex_handle(i)));
	}
	return result;
}
