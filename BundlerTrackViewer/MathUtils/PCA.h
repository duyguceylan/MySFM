/*
 *  PCA.h
 *  FunctionalModeling
 *
 *  Created by Duygu Ceylan on 9/29/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _PCA_
#define _PCA_

#include <iostream>
#include <vector>
#include "Matrix.h"

using namespace std;

template <class Type, int Dim>

class PCA
{
private:
	typedef OpenMesh::VectorT<Type,Dim> Point;
	typedef Matrix<Type,Dim,Dim> Matrix;
	
	Matrix C;					/// Covariance matrix      
	Point centroid;			// Centroid
	vector<Point> points;
	vector<Type> weights;
	
public:
	/// Constructor
	PCA();
	/// Reinitialize PCA
	void clear();
	/// Add new point
	void addPoint( Point p );
	/// Add new weighted point (don't mix with unweighted addPoint()!)
	void addPoint( Point p, Type weight);
	/// Solve eigen system
	void analyze( Point& eigenValues, Matrix& eigenVectors, Point& centroid );
};

typedef PCA<float,2> PCA2f;
typedef PCA<float,3> PCA3f;

#include "MathUtils.h"
#include "PCACode.h"

#endif
