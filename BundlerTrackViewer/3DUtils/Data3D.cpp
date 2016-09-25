/*
 *  Data3D.cpp
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/13/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "Data3D.h"
#include <iostream>

using namespace std;

int randomNumber(int maxNumberPlusOne = 0xffffffff)
{
	// rand() produces a pseudo random number between 0 and 2^15
	int ret = (rand() | 
			   ((rand() & 0x7fff) << 15 ) |
			   ((rand() & 0x7fff) << 30 )) % maxNumberPlusOne;
	return rand() % maxNumberPlusOne;
}

Data3D::Data3D()
{
	meshData.request_vertex_normals();
	meshData.request_vertex_colors();
	
	lowerCorner = Vec3f(0.0,0.0,0.0);
	upperCorner = Vec3f(0.0,0.0,0.0);
	
	OpenMesh::IO::OBJReader();
	OpenMesh::IO::OBJWriter();	
}

Data3D::Data3D(const char *filename)
{
	OpenMesh::IO::Options opt;
	bool setColor = false;
	bool setNormal = false;
	
	OpenMesh::IO::OBJReader();
	OpenMesh::IO::OBJWriter();	
	
	meshData.request_vertex_normals();
	meshData.request_vertex_colors();
	lowerCorner = Vec3f(0.0,0.0,0.0);
	upperCorner = Vec3f(0.0,0.0,0.0);
	
	string name(filename);
	modelName = name.substr(0, name.find_last_of("."));
	
	// load mesh
	if (OpenMesh::IO::read_mesh(meshData, filename, opt))
	{
		if(!opt.check(OpenMesh::IO::Options::VertexNormal))
		{
			//file does not provide vertex normals
			//compute the normals if it is a connected mesh
			if(meshData.n_faces() > 0)
			{
				meshData.request_face_normals();
				meshData.update_face_normals();
				setNormal = true;
			}
		}
		if(!opt.check(OpenMesh::IO::Options::VertexColor))
		{
			setColor = true;
		}
		
		std::cout << "mesh has " << meshData.n_vertices() << " points" << std::endl;
		
		// set bounding box
		Mesh::ConstVertexIter  vIt(meshData.vertices_begin()), 
		vEnd(meshData.vertices_end());
		
		lowerCorner = meshData.point(vIt);
		upperCorner = meshData.point(vIt);
		
		int i=0;
		for (; vIt!=vEnd; ++vIt)
		{
 			if(setNormal)
				meshData.calc_vertex_normal(vIt.handle());
			if(setColor)
				meshData.set_color(vIt, Vec3uc(175,175,175));
			
			lowerCorner.minimize(meshData.point(vIt.handle()));
			upperCorner.maximize(meshData.point(vIt.handle()));
			
			i++;
		}
		
	}
	else 
	{
		std::cerr << "[Data3D] File not found." << std::endl; 
	}
}

Data3D::~Data3D()
{
	meshData.clear();
	meshData.garbage_collection();
}

void Data3D::writeMesh(const char *filename)
{
	OpenMesh::IO::Options opt;
	opt += OpenMesh::IO::Options::VertexNormal;
	opt += OpenMesh::IO::Options::VertexColor;
	OpenMesh::IO::write_mesh(meshData, filename, opt);
}

void Data3D::setVertexPosition(int index, Vec3f pos)
{
	Mesh::VertexHandle vIt = meshData.vertex_handle(index);
	meshData.set_point(vIt, pos);
	lowerCorner.minimize(meshData.point(vIt));
	upperCorner.maximize(meshData.point(vIt));
}

void Data3D::getBoundingBox(Vec3f &lowerPt, Vec3f &upperPt)
{
	lowerPt = lowerCorner;
	upperPt	= upperCorner;
}

void Data3D::addPoint(Vec3f position, Vec3f normal, Vec3uc color)
{
	Mesh::VertexHandle vIt = meshData.new_vertex(position);
	if(meshData.has_vertex_normals())
		meshData.set_normal(vIt, normal);
	if(meshData.has_vertex_colors())
		meshData.set_color(vIt, color);
	
	if(meshData.n_vertices() == 1)
	{
		lowerCorner = position;
		upperCorner = position;
	}
	else
	{
		lowerCorner.minimize(meshData.point(vIt));
		upperCorner.maximize(meshData.point(vIt));
	}
	
	vhandles.push_back(vIt);
}

void Data3D::rigidTransformData(Matrix4f trans)
{
	Matrix3f rot = shrink4To3(trans);
	Vec3f tr(trans[3][0], trans[3][1], trans[3][2]);
	
	int noPoints = meshData.n_vertices();
	
	Vec3f pos;
	
	pos = rot * getVertexPosition(0) + tr;
	lowerCorner[0] = pos[0]; lowerCorner[1] = pos[1]; lowerCorner[2] = pos[2];
	upperCorner[0] = pos[0]; upperCorner[1] = pos[1]; upperCorner[2] = pos[2];
	
	for(int i=0; i<noPoints; i++)
	{
		pos = rot * getVertexPosition(i) + tr;
		
		lowerCorner.minimize(pos);
		upperCorner.maximize(pos);
		setVertexPosition(i, pos);
		setVertexNormal(i, rot * getVertexNormal(i));
	}
}

void Data3D::rotateData(Matrix3f rot)
{
	int noPoints = meshData.n_vertices();
	
	bool hasSearchTree = false;
	
	Vec3f pos;
	
	pos = rot * getVertexPosition(0);
	lowerCorner[0] = pos[0]; lowerCorner[1] = pos[1]; lowerCorner[2] = pos[2];
	upperCorner[0] = pos[0]; upperCorner[1] = pos[1]; upperCorner[2] = pos[2];
	
	for(int i=0; i<noPoints; i++)
	{
		pos = rot * getVertexPosition(i);
		
		lowerCorner.minimize(pos);
		upperCorner.maximize(pos);
		setVertexPosition(i, pos);
		setVertexNormal(i, rot * getVertexNormal(i));
	}
	
}

void Data3D::translateData(Vec3f trans)
{
	int dim = 3;
	int noPoints = meshData.n_vertices();
	
	Vec3f pos;
	
	pos = trans + getVertexPosition(0);
	lowerCorner[0] = pos[0]; lowerCorner[1] = pos[1]; lowerCorner[2] = pos[2];
	upperCorner[0] = pos[0]; upperCorner[1] = pos[1]; upperCorner[2] = pos[2];
	
	for(int i=0; i<noPoints; i++)
	{
		pos = getVertexPosition(i) + trans;
		lowerCorner.minimize(pos);
		upperCorner.maximize(pos);
		setVertexPosition(i, pos);
	}
	
}

void Data3D::scaleData(Vec3f scale)
{
	int noPoints = meshData.n_vertices();
	
	Vec3f pos;
	
	pos[0] = scale[0] * getVertexPosition(0)[0];
	pos[1] = scale[1] * getVertexPosition(0)[1];
	pos[2] = scale[2] * getVertexPosition(0)[2];
	lowerCorner[0] = pos[0]; lowerCorner[1] = pos[1]; lowerCorner[2] = pos[2];
	upperCorner[0] = pos[0]; upperCorner[1] = pos[1]; upperCorner[2] = pos[2];
	
	for(int i=0; i<noPoints; i++)
	{
		pos = getVertexPosition(i);
		pos[0] = scale[0] * pos[0];
		pos[1] = scale[1] * pos[1];
		pos[2] = scale[2] * pos[2];
		lowerCorner.minimize(pos);
		upperCorner.maximize(pos);
		setVertexPosition(i, pos);
	}
	
}

Data3D& Data3D::operator=(const Data3D &rhs)
{
	Vec3f pos;
	int noPoints;
	
	if(this != &rhs)
	{
		meshData = rhs.meshData;
		dataType = rhs.dataType;
		modelName = rhs.modelName;
		
		lowerCorner = rhs.lowerCorner;
		upperCorner = rhs.upperCorner;
		
	}
	
	return *this;
}

void Data3D::updateVertexNormals()
{
	Vec3f n(0.0, 0.0, 0.0);
	float norm;
	for(int i=0; i<meshData.n_vertices(); i++)
	{
		n[0]=n[1]=n[2]=0.0;
		for (Mesh::ConstVertexFaceIter vf_it=meshData.cvf_iter(meshData.vertex_handle(i)); vf_it; ++vf_it)
			n += meshData.normal(vf_it.handle());
		norm = n.length();
		if (norm != 0.0f) 
			n *= (1.0f/norm);
		meshData.set_normal(meshData.vertex_handle(i), n);
	}
}