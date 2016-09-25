/*
 *  Data3D.h
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/13/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _DATA3D_
#define _DATA3D_

#undef check

//#define id Id
#include "OpenMesh/Core/IO/MeshIO.hh"
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <OpenMesh/Core/IO/reader/OFFReader.hh>
#include <OpenMesh/Core/IO/writer/OFFWriter.hh>
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
//#undef id

#include "../MathUtils/MyMatrix.h"

using namespace std;

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

enum DATATYPE
{
	MESH = 0,
	POINTCLOUD
};

class Data3D
{
	
protected:
	
	//OpenMesh mesh
	Mesh meshData;
	
	DATATYPE dataType;
	string modelName;
	
	//bounding box
	Vec3f lowerCorner;
	Vec3f upperCorner;
	
	std::vector<Mesh::VertexHandle> vhandles;
	
public:
	Data3D(const char *filename);
	Data3D();
	~Data3D();
	
	void setType(DATATYPE t) { dataType = t;}
	
	DATATYPE getType() {return dataType; }
	string getModelName() {return modelName; }
	int getNoPoints(){return meshData.n_vertices(); }
	void getBoundingBox(Vec3f &lowerPt, Vec3f &upperPt);
	bool getHasNormal(){return meshData.has_vertex_normals(); }
	bool getHasColor(){return meshData.has_vertex_colors(); }
	Mesh* getMeshData() {return &meshData; }
	
	void writeMesh(const char *filename);
	
	void addVertexNormals() { meshData.request_vertex_normals(); }
	void addVertexColors() { meshData.request_vertex_colors(); }
	void addPoint(Vec3f position, Vec3f normal = Vec3f(0.0,0.0,0.0), Vec3uc color = Vec3uc(178,178,178));
	void updateVertexNormals();
	
	//vertex queries
	Vec3f getVertexPosition(int index) { return meshData.point(meshData.vertex_handle(index)); }
	Vec3f getVertexNormal(int index) { return meshData.normal(meshData.vertex_handle(index)); }
	Vec3uc getVertexColor(int index) { return meshData.color(meshData.vertex_handle(index));}
	
	//vertex access methods
	void setVertexPosition(int index, Vec3f pos); 
	void setVertexNormal(int index, Vec3f normal) { meshData.set_normal(meshData.vertex_handle(index), normal); }
	void setVertexColor(int index, Vec3uc col) {meshData.set_color(meshData.vertex_handle(index), col); }
	
	void rigidTransformData(Matrix4f trans);
	void rotateData(Matrix3f rot);
	void translateData(Vec3f trans);
	void scaleData(Vec3f scale);
	
	Data3D& operator=(const Data3D &rhs);
	
};

#endif