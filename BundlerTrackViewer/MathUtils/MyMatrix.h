/*
 *  Matrix.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _MATRIX_
#define	_MATRIX_

#include <iostream>

#include <OpenMesh/Core/Geometry/VectorT.hh>

typedef OpenMesh::Vec2f			Vec2f;
typedef OpenMesh::Vec2d			Vec2d;
typedef OpenMesh::Vec3f			Vec3f;
typedef OpenMesh::Vec3uc		Vec3uc;
typedef OpenMesh::Vec4f			Vec4f;
typedef OpenMesh::Vec3d			Vec3d;
typedef OpenMesh::Vec6f			Vec6f;

using namespace OpenMesh;

template <typename Type, int rows, int columns>

//this class stores data in column-major order!!!

class Matrix {
private:
	VectorT<Type, rows> theColumns[columns];
	
public:
	// ---- constructors
	Matrix<Type, rows, columns>();
	
	// ---- element access
	VectorT<Type, rows>& operator[](const unsigned &index);
	const VectorT<Type, rows>& operator[](const unsigned &index) const;
	
	//  ----  operators +, -, *, /
	inline Matrix<Type, rows, columns> operator+(const Matrix<Type, rows, columns> &op) const;
	inline Matrix<Type, rows, columns> operator-(const Matrix<Type, rows, columns> &op) const;
	inline Matrix<Type, rows, columns> operator-() const;
	Matrix<Type, rows, columns> operator*(const Type &s) const;
	
	template <int columns2>
	Matrix<Type, rows, columns2> operator*( const Matrix<Type, columns, columns2>& mat ) const;
	
	inline Matrix<Type, rows, columns> operator/(const Type &s) const;
	
	//  ---- operators +=, -=, *=, /=
	Matrix<Type, rows, columns> operator+=(const Matrix<Type, rows, columns> &op);
	inline Matrix<Type, rows, columns> operator-=(const Matrix<Type, rows, columns> &op);
	inline Matrix<Type, rows, columns> operator*=(const Type &op);
	inline Matrix<Type, rows, columns> operator*=(const Matrix<Type, rows, columns> &op);
	Matrix<Type, rows, columns> operator/=(const Type &op);
	
	//  ---- operators =, ==, !=
	inline Matrix<Type, rows, columns> operator=(const Matrix<Type, rows, columns> &op);
	inline bool operator==(const Matrix<Type, rows, columns> &op) const;
	inline bool operator!=(const Matrix<Type, rows, columns> &op) const;
	
	// ---- multiplying with vectors
	VectorT<Type, rows> operator*(const VectorT<Type, columns> &v);
	
	// ----  misc.
	inline Type* data();
	inline const Type* data() const;
	inline void loadIdentity();
	inline void changeRows(int row1, int row2);
	inline void multRow(int row, Type value);
	inline void combineRows(int row, int with, Type by);
	inline void addRows(unsigned row, unsigned with);
	Matrix<Type, columns, rows> transpose() const;
	inline void setZero();
	inline void setIdentity();
	void invertMatrix(bool *success, Type pivotEps = (Type) 1e-7f);
	Matrix<Type, rows, columns> getInverseMatrix();
	inline unsigned searchPivot(int step, bool *success, Type pivotEps);
	inline Type getDeterminant() const;
	
	/// Sort eigenvectors according to eigenvalues in descending order
	void eigenSort( VectorT<Type,rows>& eigenValues, Matrix<Type, columns, rows>& eigenVectors ) const;
	/// Compute the eigenstructure (eigenvalues and eigenvectors) of a symmetrical matrix (Numerical Recepies).
	/// The normalized eigenvectors are found as the rows of the eigenVectors-matrix.
	void computeEigenStructure( VectorT<Type,rows>& eigenValues, Matrix<Type, columns, rows>& eigenVectors ) const;
	
	/// compatibility with DynamicVector
	static inline unsigned getRowsDim() {return rows;}
	static inline unsigned getColsDim() {return columns;}
	static inline unsigned getRows()    {return rows;}
	static inline unsigned getColumns() {return columns;}
};

typedef Matrix<float, 2, 2>		Matrix2f;
typedef Matrix<float, 3, 3>		Matrix3f;
typedef Matrix<float, 4, 4>		Matrix4f;
typedef Matrix<float, 6, 6>		Matrix6f;

inline Matrix3f makeMatrix3f(float m00, float m10, float m20,
							 float m01, float m11, float m21,
							 float m02, float m12, float m22);

inline Matrix4f makeMatrix4f(float m00, float m10, float m20, float m30,
							 float m01, float m11, float m21, float m31,
							 float m02, float m12, float m22, float m32,
							 float m03, float m13, float m23, float m33);
inline Vec4f expand3To4(Vec3f v);
inline Matrix4f expand3To4(Matrix3f v);
inline Vec3f shrink4To3(Vec4f v);
inline Matrix3f shrink4To3(Matrix4f v);
inline Matrix3f makeRotX3f(double angle);
inline Matrix4f makeRotX4f(double angle);
inline Matrix3f makeRotY3f(double angle);
inline Matrix4f makeRotY4f(double angle);
inline Matrix3f makeRotZ3f(double angle);
inline Matrix4f makeRotZ4f(double angle);
inline Matrix3f makeRotVector3f(Vec3f vec, double angle);
inline Matrix4f makeRotVector4f(Vec3f vec, float angle);
inline Matrix4f makeTranslation4f(const Vec3f &t);
inline Matrix4f makeTranslation4f(float x, float y, float z);
inline Matrix3f makeScale3f(const Vec3f &s);
inline Matrix4f makeScale4f(const Vec3f &s);
inline Matrix3f makeScale3f(float x, float y, float z);
inline Matrix4f makeScale4f(float x, float y, float z);
inline Vec3f transformVector3f(Matrix4f M, const Vec3f &v);
inline Vec3f expand2To3(Vec2f v);
inline Matrix3f makeNewCoordPrjMatrix3f(Vec3f v1, Vec3f v2, Vec3f v3);
inline Matrix4f invertFrame( const Matrix4f &m );
template <class Type, int dim>
inline Matrix<Type,dim,dim> outerProduct(const VectorT<Type,dim> &v1,const VectorT<Type,dim> &v2);
inline Vec3f projectHomogeneous4To3(Vec4f v);
inline Vec2f projectHomogeneous3To2(Vec3f v);
inline Vec2f projectiveTransformVector2f(Matrix3f M, const Vec2f &v );
inline Vec3f projectiveTransformVector3f(Matrix4f M, const Vec3f &v);
inline Matrix3f calcTangentSystem( const Vec3f &normal );
inline Vec3f computeNearestPoint(const Vec3f &dir0, const Vec3f &p0, const Vec3f &dir1, const  Vec3f &p1 );
inline Vec3f computeIntersectionPoint(const Vec3f &dir0, const Vec3f &p0, const Vec3f &dir1, const Vec3f &p1 );

#include "MatrixImpl.h"

#endif