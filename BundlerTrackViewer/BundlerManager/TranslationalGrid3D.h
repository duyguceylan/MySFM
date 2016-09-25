/*
 *  TranslationalGrid3D.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 2/25/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _TRANSLATIONAL_GRID_3D_H
#define _TRANSLATIONAL_GRID_3D_H

#include "Plane.h"
#include "Mesh3D.h"
#include <vector>

using namespace std;

class TranslationalGrid3D
{
public:
	TranslationalGrid3D() { gridMesh = NULL;}
	TranslationalGrid3D(int nCols, int nRows);
	
	void setTransformations(Vec3f horTrans, Vec3f verTrans);
	void fillElementCenters(Vec3f refPoint);
	void setElement(int row, int col, vector<Vec3f> vertices);
	void setTemplateInfo(int templateImg_, int templatePlane_, Vec2f templateCenter_, Vec2i templateSize_);
	void setElementCenter(int row, int col, Vec3f pt) {elementCenters[noColumns*row+col] = pt;}
	
	void getTransformation(Vec3f &horTrans_, Vec3f &verTrans_) {horTrans_ = horizontalTrans; verTrans_ = verticalTrans; }
	int getNoElements() {return elements.size(); }
	int getNoRows() {return noRows; }
	int getNoColumns() {return noColumns; }
	vector<Vec3f> getElement(int row, int col) {/*vector<Vec3f> v; v.push_back(elementCenters[noColumns*row+col]);return v;*/ return elements[noColumns*row+col]; }
	void getTemplateInfo(int &templateImg_, int &templatePlane_, Vec2f &templateCenter_, Vec2i &templateSize_);
	Vec3f getElementCenter(int row, int col) {return elementCenters[noColumns*row+col]; }
	
	void computeGridPlaneParameters();
	void setGridPlaneParameters(Vec3f normal, float distance);
	void setGridPlane(Plane mainPlane_) {mainPlane = mainPlane_;}
	Plane getGridPlane() {return mainPlane;}
	Vec3f intersectRayWithGrid(Vec3f org, Vec3f dir);
	Vec2f transformToGridLocalCoordinates(Vec3f pos);
	Vec3f transformGridPointToGlobalCoordinateSystem(Vec2f pos);
	
	void createElementMesh();
	void createElementTexturePoints(int row, int col, int width, int height, vector<Vec3f> &texPoints);
	
	void saveGrid(const char *filename);
	void readGrid(const char *filename);
	
	void findClosestGridCell(Vec3f center, int &row, int &col);
	
	void addVisibleImage(int image) {visibleImages.push_back(image);}
	bool isVisibleInImage(int image);
	
	Plane mainPlane;
	
private:
	int noColumns;
	int noRows;
	Vec3f horizontalTrans;
	Vec3f verticalTrans;
	vector<Vec3f> elementCenters;
	vector<vector<Vec3f> > elements;
	
	//template size
	int templateImg;
	int templatePlane;
	Vec2f templateCenter;
	Vec2i templateSize;
	
	Mesh3D *gridMesh;
	
	vector<int> visibleImages;
	
};

#endif