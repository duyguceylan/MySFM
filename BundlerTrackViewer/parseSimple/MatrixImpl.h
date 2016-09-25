/*
 *  MatrixImpl.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/27/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MATRIX_IMPL_H
#define MATRIX_IMPL_H

using namespace std;

template <typename Type, int rows, int columns>
inline VectorT<Type, rows>& Matrix<Type, rows, columns>::operator[](const unsigned &index) {
	return theColumns[index];
}


template <typename Type, int rows, int columns>
inline const VectorT<Type, rows>& Matrix<Type, rows, columns>::operator[](const unsigned &index) const {
	return theColumns[index];
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator+(const Matrix<Type, rows, columns> &op) const {
	Matrix<Type, rows, columns> result;
	for (unsigned i=0; i<columns; i++) {
		result[i] = theColumns[i] + op.theColumns[i];
	}
	return result;
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator-(const Matrix<Type, rows, columns> &op) const {
	Matrix<Type, rows, columns> result;
	for (unsigned i=0; i<columns; i++) {
		result[i] = theColumns[i] - op.theColumns[i];
	}
	return result;
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator-() const {
	Matrix<Type, rows, columns> result;
	for (unsigned i=0; i<columns; i++) {
		result[i] = -theColumns[i];
	}
	return result;
}

template <typename Type, int rows, int columns>
Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator*(const Type &s) const {
	Matrix<Type, rows, columns> result;
	for (unsigned i=0; i<columns; i++) {
		result[i] = theColumns[i] * s;
	}
	return result;
}

template <typename Type, int rows, int columns>
template <int columns2>
Matrix<Type, rows, columns2> Matrix<Type, rows, columns>::operator*(const Matrix<Type, columns, columns2>& mat ) const
{
	Matrix<Type, rows, columns2> res;
	for ( int row = 0; row < rows; row++ ){
		for ( int column = 0; column < columns2; column++ ){
			res[column][row] = 0;
			for ( int k = 0; k < columns; k++ )
				res[column][row] += theColumns[k][row] * mat[column][k];
		}
	}
	return res;
}

template <typename Type, int rows, int columns>
Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator/(const Type &s) const {
	Matrix<Type, rows, columns> result;
	for (unsigned i=0; i<columns; i++) {
		result[i] = theColumns[i] / s;
	}
	return result;
}

template <typename Type, int rows, int columns>
Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator+=(const Matrix<Type, rows, columns> &op) {
	for (unsigned i=0; i<columns; i++) {
		theColumns[i] += op.theColumns[i];
	}
	return *this;
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator-=(const Matrix<Type, rows, columns> &op) {
	for (unsigned i=0; i<columns; i++) {
		theColumns[i] -= op.theColumns[i];
	}
	return *this;
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator*=(const Type &op) {
	for (unsigned i=0; i<columns; i++) {
		theColumns[i] *= op;
	}
	return *this;
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator*=(const Matrix<Type, rows, columns> &op) {
	Matrix<Type, rows, columns> result = (*this)*op;
	*this = result;
	return result;
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator/=(const Type &op) {
	for (unsigned i=0; i<columns; i++) {
		theColumns[i] /= op;
	}
	return *this;
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::operator=(const Matrix<Type, rows, columns> &op) {
	for (unsigned i=0; i<columns; i++) {
		theColumns[i] = op.theColumns[i];
	}
	return *this;
}

template <typename Type, int rows, int columns>
inline bool Matrix<Type, rows, columns>::operator==(const Matrix<Type, rows, columns> &op) const {
	for (unsigned i=0; i<columns; i++) {
		if (theColumns[i] != op.theColumns[i]) return false;
	}
	return true;
}

template <typename Type, int rows, int columns>
inline bool Matrix<Type, rows, columns>::operator!=(const Matrix<Type, rows, columns> &op) const {
	for (unsigned i=0; i<columns; i++) {
		if (theColumns[i] != op.theColumns[i]) return true;
	}
	return false;
}

template <typename Type, int rows, int columns>
VectorT<Type, rows> Matrix<Type, rows, columns>::operator*(const VectorT<Type, columns> &v) {
	VectorT<Type, rows> result;
	for (int y=0; y<rows; y++) 
	{
		result[y] = 0.0;
		for (int x=0; x<columns; x++) 
		{
			result[y] += theColumns[x][y] * v[x];
		}
	}
	return result;
}

template <typename Type, int rows, int columns>
inline Matrix<Type, rows, columns>::Matrix() {
	for (int j=0; j<columns; j++) {
		for (int i=0; i<rows; i++) {
			if (i == j)
				(theColumns[j])[i] = 1.0;
			else
				(theColumns[j])[i] = 0.0;
		}
	}
}

template <typename Type, int rows, int columns>
inline Type* Matrix<Type, rows, columns>::data() {
	return &(theColumns[0][0]);
}

template <typename Type, int rows, int columns>
inline const Type* Matrix<Type, rows, columns>::data() const {
	return &(theColumns[0][0]);
}

template <typename Type, int rows, int columns>
inline void Matrix<Type, rows, columns>::loadIdentity()
{
	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<columns; j++)
		{
			theColumns[i][j] = (Type) 0;
		}
		theColumns[i][i] = (Type) 1;
	}
}

template <typename Type, int rows, int columns>
inline void Matrix<Type, rows, columns>::changeRows(int row1, int row2) {
	VectorT<Type, columns> h;
	int i;
	for (i=0; i<columns; i++) {
		h[i] = theColumns[i][row2];
	}
	for (i=0; i<columns; i++) {
		theColumns[i][row2] = theColumns[i][row1];
	}
	for (i=0; i<columns; i++) {
		theColumns[i][row1] = h[i];
	}
}

template <typename Type, int rows, int columns>
inline void Matrix<Type, rows, columns>::multRow(int row, Type value) {
	for (int i=0; i<columns; i++) {
		theColumns[i][row] *= value;
	}
}

template <typename Type, int rows, int columns>
inline void Matrix<Type, rows, columns>::combineRows(int row, int with, Type by) {
	for (int i=0; i<columns; i++) {
		theColumns[i][row] += theColumns[i][with] * by;
	}
}

template <typename Type, int rows, int columns>
inline void Matrix<Type, rows, columns>::addRows(unsigned row, unsigned with) {
	for (unsigned i=0; i<columns; i++) {
		theColumns[i][row] += theColumns[i][with];
	}
}

template <typename Type, int rows, int columns>
Matrix<Type, columns, rows> Matrix<Type, rows, columns>::transpose() const {
	Matrix<Type, columns, rows> result;
	for (int r=0; r<rows; r++) {
		for (int c=0; c<columns; c++) {
			result[r][c] = theColumns[c][r];
		}
	}
	return result;
}

template <typename Type, int rows, int columns>
inline void Matrix<Type, rows, columns>::setZero() {
	for (int i=0; i<columns; i++) {
		for (int j=0; j<rows; j++) {
			theColumns[i][j] = (Type)0;
		}
	}
}

template <typename Type, int rows, int columns>
inline void Matrix<Type, rows, columns>::setIdentity() {
	setZero();
	unsigned minDim = rows;
	if(columns < rows)
		minDim = columns;
	for (int i=0; i<minDim; i++) {
		theColumns[i][i] = (Type)1;
	}
}

template <typename Type, int rows, int columns>
inline Type Matrix<Type, rows, columns>::getDeterminant() const {
	if ( ( rows == 2 ) && ( columns == 2 ) )
		return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
	else if ( ( rows == 3 ) && ( columns == 3 ) )
		return (*this)[0][0] * ( (*this)[1][1] * (*this)[2][2] - (*this)[1][2] * (*this)[2][1] ) -
		(*this)[0][1] * ( (*this)[1][0] * (*this)[2][2] - (*this)[1][2] * (*this)[2][0] ) +
		(*this)[0][2] * ( (*this)[1][0] * (*this)[2][1] - (*this)[1][1] * (*this)[2][0] );
	else if ( rows == 1 && columns == 1)
		return (*this)[0][0];
	else
		// limitiation: does not work for dimensions higher larger 3
		return NULL;
}

/**
 *  Helping structure for jacobi iteration
 */
#define STATIC_MATRIX_ROTATE(a,i,j,k,l) g=a[i][j];\
h=a[k][l];\
a[i][j]=g-s*(h+g*tau);\
a[k][l]=h+s*(g-h*tau);

template <typename Type, int rows, int columns>
void Matrix<Type, rows, columns>::eigenSort(VectorT<Type,rows>& eigenValues, Matrix<Type,columns,rows>& eigenVectors ) const
{
	int i,j,k;
	Type p;
	for ( i = 0; i < rows-1; i++ )
	{
		p = eigenValues[k=i];
		for( j = i+1; j < rows; j++ )
			if( eigenValues[j] >= p ) p = eigenValues[k=j];
		if( k != i ){
			eigenValues[k] = eigenValues[i];
			eigenValues[i] = p;
			for( j = 0; j < rows; j++ ){
				p = eigenVectors[j][i];
				eigenVectors[j][i] = eigenVectors[j][k];
				eigenVectors[j][k] = p;
			}
		}
	}
}

template <typename Type, int rows, int columns>
void Matrix<Type, rows, columns>::computeEigenStructure(VectorT<Type,rows>& eigenValues, Matrix<Type,columns, rows>& eigenVectors ) const
{
	if ( rows != columns ){
		cout << "ComputeEigenStructure: Only squared matrices allowed!" << endl;
		return;
	}
	
	Matrix<Type,columns, rows> a = *this;
	int j,iq,ip,i;
	Type tresh,theta,tau,t,s,h,g,c,sm = 0;
	VectorT<Type,rows> b;
	VectorT<Type,rows> z;
	for ( ip = 0; ip < rows; ip++) {
		for ( iq = 0; iq < rows; iq++ )
			eigenVectors[ip][iq] = 0.0;
		eigenVectors[ip][ip] = 1.0;
	}
	for ( ip = 0; ip < rows; ip++ ) {
		b[ip] = eigenValues[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	int nrot = 0;
	for ( i = 1; i <= 100; i++ ) {
		sm = 0;
		for ( ip = 0; ip < rows-1; ip++ ) {
			for ( iq = ip + 1; iq < rows; iq++ )
				sm += fabs( a[ip][iq] );
		}
		if ( sm == 0 ) {
			eigenSort( eigenValues, eigenVectors );
			eigenVectors = eigenVectors.transpose();
			return;
		}
		if ( i < 4 )
			tresh = 0.2f * sm / (Type)( rows * rows );
		else
			tresh = 0.0;
		for ( ip = 0; ip < rows - 1; ip++ ) {
			for ( iq = ip + 1; iq < rows; iq++ ) {
				g = (Type) 100.0f * fabs( (Type) a[ip][iq] );
				if ( ( i > 4 ) &&
					( (float)( fabs( eigenValues[ip] ) + g ) == (float)fabs( eigenValues[ip] ) ) &&
					( (float)( fabs( eigenValues[iq] ) + g ) == (float)fabs( eigenValues[iq] ) ) )
					a[ip][iq] = 0.0;
				else if ( fabs( a[ip][iq] ) > tresh ) {
					h = eigenValues[iq] - eigenValues[ip];
					if ( (float)(fabs( h ) + g ) == (float)fabs( h ) )
						t = ( a[ip][iq] ) / h;
					else {
						theta = 0.5f * h / ((Type) a[ip][iq] );
						t = 1.0f / (fabs( theta ) + sqrt( 1.0f + theta*theta ));
						if ( theta < 0.0f ) t = -t;
					}
					c = 1.0f / (Type)sqrt( 1.0f + t * t );
					s = t * c;
					tau = s / ( 1.0f + c );
					h = t * a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					eigenValues[ip] -= h;
					eigenValues[iq] += h;
					a[ip][iq] = 0.0;
					for ( j = 0; j <= ip - 1; j++ ) {
						STATIC_MATRIX_ROTATE( a, j, ip, j, iq )
					}
					for ( j = ip + 1; j <= iq - 1; j++ ) {
						STATIC_MATRIX_ROTATE( a, ip, j, j, iq )
					}
					for ( j = iq + 1; j < rows; j++ ) {
						STATIC_MATRIX_ROTATE( a, ip, j, iq, j )
					}
					for ( j = 0; j < rows; j++ ) {
						STATIC_MATRIX_ROTATE( eigenVectors, j, ip, j, iq )
					}
					++nrot;
				}
			}
		}
		for ( ip = 0; ip < rows; ip++ ) {
			b[ip] += z[ip];
			eigenValues[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	cout << "Too many iterations in routine get Eigenvalues/vectors" << endl;
}



// --- some explicit specializations that always speed up the code

template<>
inline Matrix<float, 2, 3>::Matrix() {
	theColumns[0][0] = 0.0f; theColumns[0][1] = 0.0f;
	theColumns[1][0] = 0.0f; theColumns[1][1] = 0.0f;
	theColumns[2][0] = 0.0f; theColumns[2][1] = 0.0f;
}

template<>
inline Matrix<float, 2, 2> Matrix<float, 2, 2>::operator+(const Matrix<float, 2, 2> &op) const {
	Matrix<float, 2, 2> result;
	result.theColumns[0][0] = theColumns[0][0] + op.theColumns[0][0];
	result.theColumns[0][1] = theColumns[0][1] + op.theColumns[0][1];
	result.theColumns[1][0] = theColumns[1][0] + op.theColumns[1][0];
	result.theColumns[1][1] = theColumns[1][1] + op.theColumns[1][1];
	return result;
}

template<>
inline Matrix<float, 2, 2> Matrix<float, 2, 2>::operator*(const float &s) const {
	Matrix<float, 2, 2> result;
	result.theColumns[0][0] = theColumns[0][0] * s;
	result.theColumns[0][1] = theColumns[0][1] * s;
	result.theColumns[1][0] = theColumns[1][0] * s;
	result.theColumns[1][1] = theColumns[1][1] * s;
	return result;
}

template<>
inline Matrix<float, 3, 3> Matrix<float, 3, 3>::operator+=(const Matrix<float, 3, 3> &op) {
	theColumns[0][0] += op.theColumns[0][0];
	theColumns[0][1] += op.theColumns[0][1];
	theColumns[0][2] += op.theColumns[0][2];
	
	theColumns[1][0] += op.theColumns[1][0];
	theColumns[1][1] += op.theColumns[1][1];
	theColumns[1][2] += op.theColumns[1][2];
	
	theColumns[2][0] += op.theColumns[2][0];
	theColumns[2][1] += op.theColumns[2][1];
	theColumns[2][2] += op.theColumns[2][2];
	return *this;
}

template <>
inline Matrix<float, 3, 3> Matrix<float, 3, 3>::operator*(const float &s) const {
	Matrix<float, 3, 3> result;
	result.theColumns[0][0] = theColumns[0][0] * s;
	result.theColumns[0][1] = theColumns[0][1] * s;
	result.theColumns[0][2] = theColumns[0][2] * s;
	
	result.theColumns[1][0] = theColumns[1][0] * s;
	result.theColumns[1][1] = theColumns[1][1] * s;
	result.theColumns[1][2] = theColumns[1][2] * s;
	
	result.theColumns[2][0] = theColumns[2][0] * s;
	result.theColumns[2][1] = theColumns[2][1] * s;
	result.theColumns[2][2] = theColumns[2][2] * s;
	return result;
}

template<typename Type, int rows, int columns>
void Matrix<Type, rows, columns>::invertMatrix(bool *success, Type pivotEps) {
	if(rows != columns)
	{
		cout << "invertMatrix: Matrix is not square." << endl;
		*success = false;
		return;
	}
	
	Matrix<Type, rows, columns> temp = *this;
	Matrix<Type, rows, columns> mCopy = *this; //original matrix
	
	for (int step=0; step<rows; step++) 
	{
		unsigned pivot = temp.searchPivot(step, success, pivotEps);
		if (success != NULL) {
			if (!(*success)) {
				return;
			}
		}
		if (step != pivot) {
			temp.changeRows(step, pivot);
			changeRows(step, pivot);
		}
		
		if (temp[step][step] == 0) 
		{
			cout << "invertMatrix: Matrix is singular." << endl;
			*success = false;
			*this = mCopy;
			return;
		}
		
		Type m = (Type)1.0/temp[step][step];
		temp.multRow(step, m);
		multRow(step, m);
		
		for (int row=step+1; row<rows; row++) {
			Type m = -temp[step][row];
			temp.combineRows(row, step, m);
			combineRows(row, step, m);
		}
	}
	for (int step=rows-1; step>=0; step--) 
	{
		for (int row=step-1; row>=0; row--) {
			Type m = -temp[step][row];
			//mCopy.combineRows(row, step, m);
			combineRows(row, step, m);
		}
	}
	if (success != NULL) {
		*success = true;
	}
}

template<typename Type, int rows, int columns>
inline Matrix<Type, rows, columns> Matrix<Type, rows, columns>::getInverseMatrix()
{
	Matrix<Type, rows, columns> minv;
	int indxc[rows];
	int indxr[rows];
	int ipiv[rows];
	for(int i=0;i<rows;i++)
		ipiv[i] = 0.0;
	double temp;
	
	minv = *this;
	
	for (int i = 0; i < rows; i++) {
		int irow = -1, icol = -1;
		double big = 0.;
		// Choose pivot
		for (int j = 0; j < rows; j++) {
			if (ipiv[j] != 1) {
				for (int k = 0; k < columns; k++) {
					if (ipiv[k] == 0) {
						if (fabs(minv[k][j]) >= big) {
							big = double(fabs(minv[k][j]));
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1) {
						std::cout << "ERROR: Singular matrix in MatrixInvert\n";
					}
				}
			}
		}
		++ipiv[icol];
		// Swap rows _irow_ and _icol_ for pivot
		if (irow != icol) {
			for (int k = 0; k < columns; ++k){
				temp = minv[k][irow];
				minv[k][irow] = minv[k][icol];
				minv[k][icol] = temp;
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (minv[icol][icol] == 0.){
			std::cout << "Singular matrix in MatrixInvert\n";
		}
		// Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
		double pivinv = 1.f / minv[icol][icol];
		minv[icol][icol] = 1.f;
		for (int j = 0; j < rows; j++) {
			minv[j][icol] *= pivinv;
		}
		
		// Subtract this row from others to zero out their columns
		for (int j = 0; j < rows; j++) {
			if (j != icol) {
				double save = minv[icol][j];
				minv[icol][j] = 0;
				for (int k = 0; k < columns; k++) {
					minv[k][j] -= minv[k][icol]*save;
				}
			}
		}
	}
	// Swap columns to reflect permutation
	for (int j = rows-1; j >= 0; j--) {
		if (indxr[j] != indxc[j]) {
			for (int k = 0; k < columns; k++){
				temp = minv[indxr[j]][k];
				minv[indxr[j]][k] = minv[indxc[j]][k];
				minv[indxc[j]][k] = temp;
			}
		}
	}
	return minv;
	//Matrix<Type, rows, columns> m= *this;
	//bool success;
	//m.invertMatrix(&success);
	//return m;
}

template<typename Type, int rows, int columns>
inline unsigned Matrix<Type, rows, columns>::searchPivot(int step, bool *success, Type pivotEps) {
	if(rows != columns)
	{
		cout << "searchPivot: matrix is not square." << endl;
		*success = false;
		return -1;
	}
	
	int result = step;
	Type max = fabs(theColumns[step][step]);
	for (int p = step+1; p<rows; p++) {
		if (fabs(theColumns[step][p]) > max) {
			max = fabs(theColumns[step][p]);
			result = p;
		}
	}
	if (max <= pivotEps) {
		if (success != NULL) {
			*success = false;
			return 0;
		} else {
			cout << "searchPivot: no pivot element found (matrix singular?)." << endl;
			return -1;
		}
	} else {
		if (success != NULL) {
			*success = true;
			return 0;
		}
		else {
			return -1;
		}
		
	}
}

inline Matrix3f makeMatrix3f(float m00, float m10, float m20,
							 float m01, float m11, float m21,
							 float m02, float m12, float m22)
{
	Matrix3f m;
	m[0][0] = m00; m[1][0] = m10; m[2][0] = m20;
	m[0][1] = m01; m[1][1] = m11; m[2][1] = m21;
	m[0][2] = m02; m[1][2] = m12; m[2][2] = m22;
	return m;
}

inline Matrix4f makeMatrix4f(float m00, float m10, float m20, float m30,
							 float m01, float m11, float m21, float m31,
							 float m02, float m12, float m22, float m32,
							 float m03, float m13, float m23, float m33)
{
	Matrix4f m;
	m[0][0] = m00; m[1][0] = m10; m[2][0] = m20; m[3][0] = m30;
	m[0][1] = m01; m[1][1] = m11; m[2][1] = m21; m[3][1] = m31;
	m[0][2] = m02; m[1][2] = m12; m[2][2] = m22; m[3][2] = m32;
	m[0][3] = m03; m[1][3] = m13; m[2][3] = m23; m[3][3] = m33;
	return m;
}

inline Matrix4f makeMatrix4f(Vec4f c1, Vec4f c2, Vec4f c3, Vec4f c4)
{
	Matrix4f m;
	m[0] = c1;
	m[1] = c2;
	m[2] = c3;
	m[3] = c4;
	return m;
}

inline Vec4f expand3To4(Vec3f v)
{
	return Vec4f(v[0],v[1],v[2],1.0);
}

inline Matrix4f expand3To4(Matrix3f v) {
	Matrix4f result;
	for (int x = 0; x<3; x++) {
		for (int y = 0; y<3; y++) {
			result[x][y] = v[x][y];
		}
	}
	for (int i = 0; i<3; i++) {
		result[3][i] = 0.0;
	}
	for (int i = 0; i<3; i++) {
		result[i][3] = 0.0;
	}
	result[3][3] = 1.0;
	return result;
}

inline Vec3f shrink4To3(Vec4f v) {
	return Vec3f(v[0], v[1], v[2]);
}

inline Matrix3f shrink4To3(Matrix4f v) {
	Matrix3f result;
	for (int x = 0; x<3; x++) {
		for (int y = 0; y<3; y++) {
			result[x][y] = v[x][y];
		}
	}
	return result;
}

inline Vec2f shrink3To2(Vec3f v)
{
	return Vec2f(v[0], v[1]);
}

//angle is in radians
inline Matrix3f makeRotX3f(double angle) {
	return makeMatrix3f (
						 1.0,        0.0,        0.0,
						 
						 0.0, cos(angle), sin(angle),
						 
						 0.0,-sin(angle), cos(angle)
						 );
}

//angle is in radians
inline Matrix4f makeRotX4f(double angle) {
	return expand3To4(makeRotX3f(angle));
}

//angle is in radians
inline Matrix3f makeRotY3f(double angle) {
	return makeMatrix3f (
						 cos(angle),        0.0, sin(angle),
						 
						 0.0,        1.0,        0.0,
						 
						 -sin(angle),        0.0, cos(angle)
						 );
}

//angle is in radians
inline Matrix4f makeRotY4f(double angle) {
	return expand3To4(makeRotY3f(angle));
}

//angle is in radians
inline Matrix3f makeRotZ3f(double angle) {
	return makeMatrix3f (
						 cos(angle), sin(angle),        0.0,
						 
						 -sin(angle), cos(angle),        0.0,
						 
						 0.0,        0.0,        1.0
						 );
}

//angle is in radians
inline Matrix4f makeRotZ4f(double angle) {
	return expand3To4(makeRotZ3f(angle));
}

//angle is in radians
inline Matrix3f makeRotVector3f(Vec3f vec, double angle) {
	float cosA    = cos(angle),
	sinA    = sin(angle),
	invCosA = 1.0f - cosA;
	return makeMatrix3f( invCosA*vec[0]*vec[0] + cosA,
                        invCosA*vec[0]*vec[1] + sinA*vec[2],
                        invCosA*vec[0]*vec[2] - sinA*vec[1],
						
                        invCosA*vec[0]*vec[1] - sinA*vec[2],
                        invCosA*vec[1]*vec[1] + cosA,
                        invCosA*vec[1]*vec[2] + sinA*vec[0],
						
                        invCosA*vec[0]*vec[2] + sinA*vec[1],
                        invCosA*vec[1]*vec[2] - sinA*vec[0],
                        invCosA*vec[2]*vec[2] + cosA).transpose();
}

inline Matrix3f makeRotEulerAngles(float alpha, float beta, float gamma)
{
	Matrix3f R;
	
	float sa = sin(alpha);
	float ca = sqrt(1-sa*sa);
	float sb = sin(beta);
	float cb = sqrt(1-sb*sb);
	float sr = sin(gamma);
	float cr = sqrt(1-sr*sr);
	
	R[0][0] = cb*cr;
	R[1][0] = -cb*sr;
	R[2][0] = sb;
	
	R[0][1] = sa*sb*cr + ca*sr;
	R[1][1] = -sa*sb*sr + ca*cr;
	R[2][1] = -sa*cb;
	
	R[0][2] = -ca*sb*cr + sa*sr;
	R[1][2] = ca*sb*sr + sa*cr;
	R[2][2] = ca*cb;
	
	return R;
}


inline Matrix4f makeRotVector4f(Vec3f vec, float angle) {
	return expand3To4(makeRotVector3f(vec, angle));
}

inline Matrix4f makeTranslation4f(const Vec3f &t) {
	return makeMatrix4f (
						 1.0,  0.0,  0.0, t[0],
						 0.0,  1.0,  0.0, t[1],
						 0.0,  0.0,  1.0, t[2],
						 0.0,  0.0,  0.0,  1.0
						 );
}

inline Matrix4f makeTranslation4f(float x, float y, float z) {
	return makeMatrix4f (
						 1.0,  0.0,  0.0,   x,
						 0.0,  1.0,  0.0,   y,
						 0.0,  0.0,  1.0,   z,
						 0.0,  0.0,  0.0,  1.0
						 );
}

inline Matrix3f makeScale3f(const Vec3f &s) {
	return makeMatrix3f (
						 s[0],  0.0,  0.0,
						 0.0, s[1],  0.0,
						 0.0,  0.0, s[2]
						 );
}

inline Matrix4f makeScale4f(const Vec3f &s) {
	return expand3To4(makeScale3f(s));
}

inline Matrix3f makeScale3f(float x, float y, float z) {
	return makeScale3f(Vec3f(x,y,z));
}


inline Vec3f transformVector3f(Matrix4f M, const Vec3f &v) {
	Vec4f v4 = expand3To4(v);
	Vec4f trans = M * v4;
	return shrink4To3(trans);
}

inline Vec3f expand2To3(Vec2f v) {
	return Vec3f(v[0], v[1],  1.0f);
}

inline Matrix3f makeNewCoordPrjMatrix3f(Vec3f v1, Vec3f v2, Vec3f v3) {
	Matrix3f result;
	result[0] = v1;
	result[1] = v2;
	result[2] = v3;
	return result;
}

inline Matrix4f invertFrame( const Matrix4f &m )
{
	Matrix4f A = expand3To4( shrink4To3(m).transpose());
	Matrix4f B = makeTranslation4f( -shrink4To3(m[3]));
	return A*B;
}

template <class Type, int dim>
inline Matrix<Type,dim,dim> outerProduct(const VectorT<Type,dim> &v1,const VectorT<Type,dim> &v2)
{
	Matrix<Type,dim,dim> ret;
	for( int i=0;i<dim;i++)
	{
		for( int j=0;j<dim;j++ )
		{
			ret[i][j] = v1[i]*v2[j];
		}
	}
	return ret;
}

inline Vec3f projectHomogeneous4To3(Vec4f v) {
	return Vec3f(v[0]/v[3], v[1]/v[3], v[2]/v[3]);
}

inline Vec2f projectHomogeneous3To2(Vec3f v) {
	return Vec2f(v[0]/v[2], v[1]/v[2]);
}

inline Vec2f projectiveTransformVector2f(Matrix3f M, const Vec2f &v )
{
	return projectHomogeneous3To2(M*(expand2To3(v)));
}

inline Vec3f projectiveTransformVector3f(Matrix4f M, const Vec3f &v) {
	
	return projectHomogeneous4To3(M*(expand3To4(v)));
}

inline Matrix3f calcTangentSystem( const Vec3f &normal )
{
	Vec3f n = normal;
	n.normalize();
	Vec3f u,v;
	
	if( fabs(n[0]) > 0.5f )
		u = cross(n, Vec3f(0,1,0));
	else
		u = cross(n, Vec3f(1,0,0));
	v= cross(u,n);
	u = cross(v,n);
	u.normalize();
	v.normalize();
	return makeMatrix3f( u[0], v[0], n[0], 
						u[1], v[1], n[1],
						u[2], v[2], n[2] );
}

inline Vec3f computeNearestPoint(const Vec3f &dir0, const Vec3f &p0, const Vec3f &dir1, const  Vec3f &p1 )
{
	// transform line 0 into local frame of line 1
	Matrix3f frame1 = calcTangentSystem(dir1);
	Vec3f g = frame1.transpose() * dir0;
	Vec3f g0 = frame1.transpose() * (p0 - p1);
	g[2] = 0;
	g0[2] = 0;
	
	float t = dot(g, g0) / dot(g, g);
	
	return p0 - dir0*t;
}

inline Vec3f computeIntersectionPoint(const Vec3f &dir0, const Vec3f &p0, const Vec3f &dir1, const Vec3f &p1 )
{
	return computeNearestPoint(dir0, p0, dir1, p1)*0.5f + computeNearestPoint(dir1, p1, dir0, p0) * 0.5f;
}


#endif