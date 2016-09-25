
/*
 *  Mesh3D.cpp
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/13/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <OpenGL/OpenGL.h>
#include "Mesh3D.h"
#include <fstream>

Mesh3D::Mesh3D() : Data3D() 
{
	dataType = MESH;
}

Mesh3D::Mesh3D(const char *filename) : Data3D(filename) 
{
	dataType = MESH;
}

int Mesh3D::getVertexOfFace(int faceId, int index) 
{
	Mesh::FaceHandle fh = meshData.face_handle(faceId);
	Mesh::FaceVertexIter fv_it=meshData.fv_iter(fh);
	int v1 = fv_it.handle().idx(); ++fv_it;
	int v2 = fv_it.handle().idx(); ++fv_it;
	int v3 = fv_it.handle().idx();
	if(index == 0)
		return v1;
	else if(index == 1)
		return v2;
	else if(index == 2)
		return v3;
}

void Mesh3D::addFace(int vId1, int vId2, int vId3)
{
	std::vector<Mesh::VertexHandle>  face_vhandles;
	face_vhandles.push_back(meshData.vertex_handle(vId1));
	face_vhandles.push_back(meshData.vertex_handle(vId2));
	face_vhandles.push_back(meshData.vertex_handle(vId3));
	meshData.add_face(face_vhandles);
}

void Mesh3D::addFace(int vId1, int vId2, int vId3, int vId4)
{
	std::vector<Mesh::VertexHandle>  face_vhandles;
	face_vhandles.push_back(meshData.vertex_handle(vId1));
	face_vhandles.push_back(meshData.vertex_handle(vId2));
	face_vhandles.push_back(meshData.vertex_handle(vId3));
	face_vhandles.push_back(meshData.vertex_handle(vId4));
	
	meshData.add_face(face_vhandles);
	
}

Vec3f Mesh3D::samplePointOnFace(int index)
{
	double	rnd_1, rnd_2;
	
	// calculate two edges of face
	Mesh::FaceHandle fHandle = meshData.face_handle(index);
	Mesh::HalfedgeHandle hh = meshData.halfedge_handle(fHandle);
	Vec3f p1 = meshData.point(meshData.to_vertex_handle(hh));
	hh = meshData.next_halfedge_handle(hh);
	Vec3f p2 = meshData.point(meshData.to_vertex_handle(hh));
	hh = meshData.next_halfedge_handle(hh);
	Vec3f p3 = meshData.point(meshData.to_vertex_handle(hh));
	
	Vec3f v1 = p2 - p1;	// v1 = p2-p1 
	Vec3f v2 = p3 - p1;	// v2 = p3-p1
	
	// choose two random numbers.
	rnd_1 = (double)rand() / (double)RAND_MAX;
	rnd_2 = (double)rand() / (double)RAND_MAX;
	
	if(rnd_1 + rnd_2 > 1.0)
	{
		rnd_1 = 1.0 - rnd_1;
		rnd_2 = 1.0 - rnd_2;
	}
	
	Vec3f tmp = (v1 * rnd_1) + (v2 * rnd_2);
	
	return Vec3f(tmp[0] + p1[0], tmp[1] + p1[1], tmp[2] + p1[2]);
}

Vec3f Mesh3D::computeBarycenterCoord(int index, Vec3f pos)
{
	double l0, l1, l2, A, B, C, D, E, F, G, H, I;
	Mesh::FaceHandle fHandle = meshData.face_handle(index);
	Mesh::HalfedgeHandle hh = meshData.halfedge_handle(fHandle);
	Vec3f p1 = meshData.point(meshData.to_vertex_handle(hh));
	hh = meshData.next_halfedge_handle(hh);
	Vec3f p2 = meshData.point(meshData.to_vertex_handle(hh));
	hh = meshData.next_halfedge_handle(hh);
	Vec3f p3 = meshData.point(meshData.to_vertex_handle(hh));
	
	A = p1[0] - p3[0];
	B = p2[0] - p3[0];
	C = p3[0] - pos[0];
	D = p1[1] - p3[1];
	E = p2[1] - p3[1];
	F = p3[1] - pos[1];
	G = p1[2] - p3[2];
	H = p2[2] - p3[2];
	I = p3[2] - pos[2];
	l0 = (B*(F+I) - C*(E+H))/(A*(E+H) - B*(D+G));
	l1 = (A*(F+I) - C*(D+G))/(B*(D+G) - A*(E+H));
	l2 = 1.0 - l0 - l1;
	
	return Vec3f(l0, l1, l2);
}

Vec3f Mesh3D::computeFromBarycenterCoord(int index, Vec3f pos)
{
	double x0,y0,z0;
	Mesh::FaceHandle fHandle = meshData.face_handle(index);
	Mesh::HalfedgeHandle hh = meshData.halfedge_handle(fHandle);
	Vec3f p1 = meshData.point(meshData.to_vertex_handle(hh));
	hh = meshData.next_halfedge_handle(hh);
	Vec3f p2 = meshData.point(meshData.to_vertex_handle(hh));
	hh = meshData.next_halfedge_handle(hh);
	Vec3f p3 = meshData.point(meshData.to_vertex_handle(hh));
	
	x0 = pos[0]*p1[0] + pos[1]*p2[0] + pos[2]*p3[0];
	y0 = pos[0]*p1[1] + pos[1]*p2[1] + pos[2]*p3[1];
	z0 = pos[0]*p1[2] + pos[1]*p2[2] + pos[2]*p3[2];
	return Vec3f(x0, y0, z0);
}

void Mesh3D::computeKRingNeighborhood(int k, int index, vector<int> &indices)
{
	Mesh::VertexHandle vh;
	Mesh::VertexVertexIter vvIt;
	
	vh = meshData.vertex_handle(index);
	vvIt = meshData.vv_iter(vh);
	
	if(k==1)
	{
		for(;vvIt;++vvIt)
		{
			indices.push_back(vvIt.handle().idx());
		}
	}
	else
	{
		indices.push_back(index);
		for(;vvIt;++vvIt)
		{
			computeKRingNeighborhood(k-1, vvIt.handle().idx(), indices);
		}
	}
}