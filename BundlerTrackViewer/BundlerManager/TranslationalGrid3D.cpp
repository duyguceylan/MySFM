/*
 *  TranslationalGrid3D.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 2/25/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include "TranslationalGrid3D.h"
#include "../Triangulation/TriangulateEC.h"

TranslationalGrid3D::TranslationalGrid3D(int nCols, int nRows)
{
	noColumns = nCols;
	noRows = nRows;
	elements.clear();
	elements.resize(noColumns*noRows);
	gridMesh = NULL;
}

void TranslationalGrid3D::setTransformations(Vec3f horTrans, Vec3f verTrans)
{
	horizontalTrans = horTrans;
	verticalTrans = verTrans;
}

void TranslationalGrid3D::setTemplateInfo(int templateImg_, int templatePlane_, Vec2f templateCenter_, Vec2i templateSize_)
{
	templateImg = templateImg_;
	templatePlane = templatePlane_;
	templateCenter = templateCenter_;
	templateSize = templateSize_;
}

void TranslationalGrid3D::getTemplateInfo(int &templateImg_, int &templatePlane_, Vec2f &templateCenter_, Vec2i &templateSize_)
{
	templateImg_ = templateImg;
	templatePlane_ = templatePlane;
	templateCenter_ = templateCenter;
	templateSize_ = templateSize;
}

void TranslationalGrid3D::fillElementCenters(Vec3f refPoint)
{
	elementCenters.clear();
	for(int i=0; i<noRows; i++)
	{
		for(int j=0; j<noColumns; j++)
		{
			elementCenters.push_back(refPoint + horizontalTrans*j + verticalTrans*i);
		}
	}
}

void TranslationalGrid3D::saveGrid(const char *filename)
{
	ofstream fout(filename, ios::out);
	fout << noColumns << " " << noRows << endl;
	fout << horizontalTrans[0] << " " << horizontalTrans[1] << " " << horizontalTrans[2] << " "
		 << verticalTrans[0] << " " << verticalTrans[1] << " " << verticalTrans[2] << endl;
	for(int i=0; i<elementCenters.size(); i++)
	{
		fout << elementCenters[i][0] << " " << elementCenters[i][1] << " " << elementCenters[i][2] << endl;
	}
	fout.close();
}

void TranslationalGrid3D::readGrid(const char *filename)
{
	ifstream fin(filename, ios::in);
	if(fin)
	{
		fin >> noColumns >> noRows;
		fin >> horizontalTrans[0] >> horizontalTrans[1] >> horizontalTrans[2];
		fin >> verticalTrans[0] >> verticalTrans[1] >> verticalTrans[2];
		for(int i=0; i<noColumns*noRows; i++)
		{
			Vec3f center;
			fin >> center[0] >> center[1] >> center[2];
			elementCenters.push_back(center);
		}
		elements.clear();
		elements.resize(noColumns*noRows);
	}
	else
	{
		printf("Grid file not found.\n");
	}
}

void TranslationalGrid3D::computeGridPlaneParameters()
{
	Vec3f xDir = horizontalTrans;
	Vec3f yDir = verticalTrans;
	
	if(elementCenters.size() > 0)
	{
		Vec3f planeNormal = (cross(xDir.normalize(), yDir.normalize())).normalize();
		planeNormal = planeNormal.normalize();
		
		float planeOffset = 0.0 - dot(planeNormal, elementCenters[0]);
		
		mainPlane = Plane(planeNormal, planeOffset);
		mainPlane.setAxis(xDir, yDir);
	}
	else
	{
		printf("No element on grid specified.\n");
	}
}

void TranslationalGrid3D::setGridPlaneParameters(Vec3f normal, float distance)
{
	mainPlane = Plane(normal, distance);
}

Vec3f TranslationalGrid3D::intersectRayWithGrid(Vec3f org, Vec3f dir)
{
	Vec3f intersection;
	mainPlane.findPlaneRayIntersection(org, dir, intersection);
	return intersection;
}

Vec2f TranslationalGrid3D::transformToGridLocalCoordinates(Vec3f pos)
{
	return mainPlane.transformPointToLocalCoordinateSystem(pos);
}

Vec3f TranslationalGrid3D::transformGridPointToGlobalCoordinateSystem(Vec2f pos)
{
	return mainPlane.transformPointToGlobalCoordinateSystem(pos);
}

void TranslationalGrid3D::setElement(int row, int column, vector<Vec3f> vertices)
{
	vector<Vec3f> tempVertices;
	for(int i=0; i<vertices.size(); i++)
	{
		tempVertices.push_back(vertices[i]-verticalTrans*row-horizontalTrans*column);
	}
	
	
	for(int i=0; i<noRows; i++)
	{
		for(int j=0; j<noColumns; j++)
		{
			for(int k=0; k<vertices.size(); k++)
				elements[i*noColumns+j].push_back(tempVertices[k]+verticalTrans*i+horizontalTrans*j);
		}
	}
}

void TranslationalGrid3D::createElementMesh()
{
	if(gridMesh == NULL)
	{
		gridMesh = new Mesh3D();
		//gridMesh->getMeshData()->request_vertex_status();
		//gridMesh->getMeshData()->request_edge_status();
		//gridMesh->getMeshData()->request_face_status();
	}
	
	int previousNoPoints;
	float epsilon = 1E-5;
	
	for(int i=0; i<noRows; i++)
	{
		for(int j=0; j<noColumns; j++)
		{
			std::vector<int> outerIndices;
			vector<Vec2f> vertices;
			std::vector<int> mTriangles;
			previousNoPoints = gridMesh->getNoPoints();
			printf("******\n");
			for(int v=0; v<elements[i*noColumns+j].size(); v++)
			{
				gridMesh->addPoint(elements[i*noColumns+j][v]);
				vertices.push_back(transformToGridLocalCoordinates(elements[i*noColumns+j][v]));
				outerIndices.push_back(v);
				printf("element:%f %f %f\n", elements[i*noColumns+j][v][0], elements[i*noColumns+j][v][1], elements[i*noColumns+j][v][2]);
				printf("planar:%f %f\n", vertices[v][0], vertices[v][1]);
			}
			
			gridMesh->addFace(previousNoPoints, previousNoPoints+1, previousNoPoints+2);
			gridMesh->addFace(previousNoPoints, previousNoPoints+2, previousNoPoints+3);
			
			//TriangulateEC(vertices, Query::QT_REAL, epsilon, outerIndices, mTriangles);
			//for(int k=0; k<mTriangles.size(); k+=3)
			//{
			//	gridMesh->addFace(mTriangles[k]+previousNoPoints, mTriangles[k+1]+previousNoPoints, mTriangles[k+2]+previousNoPoints);
			//}
		}
	}
	
	std::vector<int> outerIndices;
	vector<Vec2f> vertices;
	std::vector<int> mTriangles;
	for(int j=0; j<noRows; j++)
	{
		previousNoPoints = gridMesh->getNoPoints();
		for(int k=0; k<elements[j*noColumns].size(); k++)
		{
			gridMesh->addPoint(elements[j*noColumns][k]-horizontalTrans);
			vertices.push_back(transformToGridLocalCoordinates(elements[j*noColumns][k]-horizontalTrans));
			outerIndices.push_back(k);
		}
		
		gridMesh->addFace(previousNoPoints, previousNoPoints+1, previousNoPoints+2);
		gridMesh->addFace(previousNoPoints, previousNoPoints+2, previousNoPoints+3);
		//TriangulateEC(vertices, Query::QT_REAL, epsilon, outerIndices, mTriangles);
		//for(int k=0; k<mTriangles.size(); k+=3)
		//{
		//	gridMesh->addFace(mTriangles[k]+previousNoPoints, mTriangles[k+1]+previousNoPoints, mTriangles[k+2]+previousNoPoints);
		//}
	}
	
	
	for(int j=0; j<1; j++)
	{
		previousNoPoints = gridMesh->getNoPoints();
		for(int k=0; k<elements[j*noColumns].size(); k++)
		{
			gridMesh->addPoint(elements[j*noColumns][k]-horizontalTrans-verticalTrans);
			vertices.push_back(transformToGridLocalCoordinates(elements[j*noColumns][k]-horizontalTrans-verticalTrans));
			outerIndices.push_back(k);
		}
		
		gridMesh->addFace(previousNoPoints, previousNoPoints+1, previousNoPoints+2);
		gridMesh->addFace(previousNoPoints, previousNoPoints+2, previousNoPoints+3);
		
		previousNoPoints = gridMesh->getNoPoints();
		for(int k=0; k<elements[j*noColumns].size(); k++)
		{
			gridMesh->addPoint(elements[j*noColumns][k]-verticalTrans);
			vertices.push_back(transformToGridLocalCoordinates(elements[j*noColumns][k]-verticalTrans));
			outerIndices.push_back(k);
		}
		
		gridMesh->addFace(previousNoPoints, previousNoPoints+1, previousNoPoints+2);
		gridMesh->addFace(previousNoPoints, previousNoPoints+2, previousNoPoints+3);
		
		previousNoPoints = gridMesh->getNoPoints();
		for(int k=0; k<elements[j*noColumns].size(); k++)
		{
			gridMesh->addPoint(elements[j*noColumns][k]-verticalTrans+horizontalTrans);
			vertices.push_back(transformToGridLocalCoordinates(elements[j*noColumns][k]-verticalTrans+horizontalTrans));
			outerIndices.push_back(k);
		}
		
		gridMesh->addFace(previousNoPoints, previousNoPoints+1, previousNoPoints+2);
		gridMesh->addFace(previousNoPoints, previousNoPoints+2, previousNoPoints+3);
		
		previousNoPoints = gridMesh->getNoPoints();
		for(int k=0; k<elements[j*noColumns].size(); k++)
		{
			gridMesh->addPoint(elements[j*noColumns][k]-verticalTrans+horizontalTrans*2);
			vertices.push_back(transformToGridLocalCoordinates(elements[j*noColumns][k]-verticalTrans+horizontalTrans*2));
			outerIndices.push_back(k);
		}
		
		gridMesh->addFace(previousNoPoints, previousNoPoints+1, previousNoPoints+2);
		gridMesh->addFace(previousNoPoints, previousNoPoints+2, previousNoPoints+3);
	}
	
	gridMesh->writeMesh("/Users/ceylan/Desktop/planes.obj");
	
}

void TranslationalGrid3D::createElementTexturePoints(int row, int col, int width, int height, vector<Vec3f> &texPoints)
{
	Vec2f minVertex, maxVertex;
	for(int v=0; v<elements[row*noColumns+col].size(); v++)
	{
		Vec2f pos = transformToGridLocalCoordinates(elements[row*noColumns+col][v]);
		if(v==0)
		{
			minVertex = pos;
			maxVertex = pos;
		}
		minVertex.minimize(pos);
		maxVertex.maximize(pos);
	}
	
	Vec2f xDir(maxVertex[0]-minVertex[0], 0.0);
	Vec2f yDir(0.0, maxVertex[1]-minVertex[1]);
	xDir /= width;
	yDir /= height;
	
	int patchSize = 3;
	for(int y=0; y<height; y++)
	{
		for(int x=0; x<width; x++)
		{
			Vec2f pos = minVertex + xDir*x + yDir*y;
			texPoints.push_back(transformGridPointToGlobalCoordinateSystem(pos)); 
		}
	}
}

void TranslationalGrid3D::findClosestGridCell(Vec3f center, int &row, int &col)
{
	float minDist;
	int minIndex =-1;
	for(int r=0; r<noRows; r++)
	{
		for(int c=0; c<noColumns; c++)
		{
			if(minIndex == -1)
			{
				minIndex = r*noColumns + c;
				row = r;
				col = c;
				minDist = (elementCenters[minIndex] - center).length();
			}
			else if((elementCenters[r*noColumns + c] - center).length() < minDist)
			{
				minDist = (elementCenters[r*noColumns + c] - center).length();
				row = r;
				col = c;
			}
		}
	}
}

bool TranslationalGrid3D::isVisibleInImage(int image)
{
	if(find(visibleImages.begin(), visibleImages.end(), image) == visibleImages.end())
		return false;
	else
		return true;
}