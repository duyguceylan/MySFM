/*
 *  Mesh3D.h
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/13/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _MESH3D_
#define _MESH3D_

#include "Data3D.h"

class Mesh3D : public Data3D
{
private:
	
public:
	Mesh3D();
	Mesh3D(const char *filename);
	
	void addFaceNormals() { meshData.request_face_normals(); }
	void addFace(int vId1, int vId2, int vId3);
	void addFace(int vId1, int vId2, int vId3, int vId4);
	//face queries
	int getNoFaces() {return meshData.n_faces(); }
	int getVertexOfFace(int faceId, int index);
	void updateFaceNormals() { meshData.request_face_normals(); meshData.update_face_normals(); }
	
	Vec3f samplePointOnFace(int index);
	//compute barycenter coordinates for the euclidean point
	Vec3f computeBarycenterCoord(int index, Vec3f pos);
	//compute euclidean coordinates for the barycenter point
	Vec3f computeFromBarycenterCoord(int index, Vec3f);
	void computeKRingNeighborhood(int k, int index, vector<int> &indices);
};

#endif