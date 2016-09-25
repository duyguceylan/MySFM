/*
 *  Camera.cpp
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <fstream>

#include "PerspectiveCamera.h"
#include "../MathUtils/MathUtils.h"

#include <cstdlib>

using namespace std;

PerspectiveCamera::PerspectiveCamera(void) {
	axesScale = 1.0f;
	maxLevel = 1;
	// initialize camera
	intrinsics.resize(6);
	extrinsics.resize(6);
	
	projection.resize(maxLevel);
	for (int level = 0; level < maxLevel; ++level)
	{
		Matrix4f m;
		m.loadIdentity();
		projection[level] = m;
	}
}

PerspectiveCamera::~PerspectiveCamera() {}

void PerspectiveCamera::init(const std::string cname_, const int maxLevel_) 
{
	camParamFilename = cname_;
	maxLevel = maxLevel_;
	
	// initialize camera
	intrinsics.resize(6);
	extrinsics.resize(6);
	
	ifstream ifstr;
	ifstr.open(cname_.c_str());
	
	string header;
	ifstr >> header;
	if (header != "CONTOUR")
	{
		cerr << "Unrecognizable txt format" << endl;
		exit (1);
	}
	
	for (int i = 0; i < 6; ++i)
		ifstr >> intrinsics[i];
	for (int i = 0; i < 6; ++i)
		ifstr >> extrinsics[i];
	
	ifstr.close();
	
	//----------------------------------------------------------------------
	projection.resize(maxLevel);
	for (int level = 0; level < maxLevel; ++level)
	{
		Matrix4f m;
		m.loadIdentity();
		projection[level] = m;
	}
	
	updateCamera();
}

void PerspectiveCamera::setProjectionMatrix(Matrix4f mat)
{
	projection[0] = mat;
	for (int level = 1; level < maxLevel; ++level) 
	{
		projection[level] = projection[level - 1];
		
		//halve the first and second rows
		for(int i=0; i<4; i++)
		{
			projection[level][i][0] /= 2.0;
			projection[level][i][1] /= 2.0;
		}
	}
	
	//----------------------------------------------------------------------
	//oaxis is the third row
	oaxis[0] = (projection[0])[0][2];
	oaxis[1] = (projection[0])[1][2];
	oaxis[2] = (projection[0])[2][2];
	oaxis[3] = 0.0;
	
	const float ftmp = oaxis.norm();
	oaxis[3] = (projection[0])[3][2];
	oaxis /= ftmp;
	
	center = getOpticalCenter();
	
	zaxis = Vec3f(oaxis[0], oaxis[1], oaxis[2]);
	xaxis = Vec3f((projection[0])[0][0],
				  (projection[0])[1][0],
				  (projection[0])[2][0]);
	yaxis = cross(zaxis, xaxis);
	yaxis.normalize();
	xaxis = cross(yaxis, zaxis);
	
	Vec4f xaxis_tmp, yaxis_tmp;
	xaxis_tmp[0] = (projection[0])[0][0]; xaxis_tmp[1] = (projection[0])[1][0]; xaxis_tmp[2] = (projection[0])[2][0];
	xaxis_tmp[3] = 0.0f; 
	yaxis_tmp[0] = (projection[0])[0][1]; yaxis_tmp[1] = (projection[0])[1][1]; yaxis_tmp[2] = (projection[0])[2][1];
	yaxis_tmp[3] = 0.0f;  
	//float ftmp2 = (xaxis_tmp.norm() + yaxis_tmp.norm()) / 2.0f;
	float ftmp2 = (dot(xaxis_tmp, expand3To4(xaxis)) + dot(yaxis_tmp, expand3To4(yaxis))) / 2.0;
	if (ftmp2 == 0.0f)
		ftmp2 = 1.0f;
	ipscale = ftmp2;
	
}

void PerspectiveCamera::updateProjection(void) 
{
	// Set bottom level
	setProjection(intrinsics, extrinsics, projection[0]);
	
	for (int level = 1; level < maxLevel; ++level) 
	{
		projection[level] = projection[level - 1];
		
		//halve the first and second rows
		for(int i=0; i<4; i++)
		{
			projection[level][i][0] /= 2.0;
			projection[level][i][1] /= 2.0;
		}
	}
}

void PerspectiveCamera::write(const std::string file) 
{
	ofstream ofstr;
	ofstr.open(file.c_str());
	ofstr << "CONTOUR" << endl
	<< intrinsics[0] << ' ' << intrinsics[1] << ' '
	<< intrinsics[2] << ' ' << intrinsics[3] << endl
	<< intrinsics[4] << ' ' << intrinsics[5] << ' '
	<< extrinsics[0] << ' ' << extrinsics[1] << endl
	<< extrinsics[2] << ' ' << extrinsics[3] << ' '
	<< extrinsics[4] << ' ' << extrinsics[5] << endl;
	
	ofstr.close();
}

void PerspectiveCamera::updateCamera(void) 
{
	updateProjection();
	
	//----------------------------------------------------------------------
	//oaxis is the third row
	oaxis[0] = (projection[0])[0][2];
	oaxis[1] = (projection[0])[1][2];
	oaxis[2] = (projection[0])[2][2];
	oaxis[3] = 0.0;
	
	const float ftmp = oaxis.norm();
	oaxis[3] = (projection[0])[3][2];
	oaxis /= ftmp;
	
	center = getOpticalCenter();
	
	zaxis = Vec3f(oaxis[0], oaxis[1], oaxis[2]);
	xaxis = Vec3f((projection[0])[0][0],
				  (projection[0])[1][0],
				  (projection[0])[2][0]);
	yaxis = cross(zaxis, xaxis);
	yaxis.normalize();
	xaxis = cross(yaxis, zaxis);
	
	Vec4f xaxis_tmp, yaxis_tmp;
	xaxis_tmp[0] = (projection[0])[0][0]; xaxis_tmp[1] = (projection[0])[1][0]; xaxis_tmp[2] = (projection[0])[2][0];
	xaxis_tmp[3] = 0.0f; 
	yaxis_tmp[0] = (projection[0])[0][1]; yaxis_tmp[1] = (projection[0])[1][1]; yaxis_tmp[2] = (projection[0])[2][1];
	yaxis_tmp[3] = 0.0f;  
	//float ftmp2 = (xaxis_tmp.norm() + yaxis_tmp.norm()) / 2.0f;
	float ftmp2 = (dot(xaxis_tmp, expand3To4(xaxis)) + dot(yaxis_tmp, expand3To4(yaxis))) / 2.0;
	if (ftmp2 == 0.0f)
		ftmp2 = 1.0f;
	ipscale = ftmp2;
}

Vec4f PerspectiveCamera::getOpticalCenter(void) const 
{
	// orthographic case
	Vec4f ans;
	if ((projection[0])[0][2] == 0.0 && (projection[0])[1][2] == 0.0 &&
		(projection[0])[2][2] == 0.0) 
	{
		Vec3f vtmp[2];
		for (int i = 0; i < 2; ++i)
			for (int y = 0; y < 3; ++y)
				vtmp[i][y] = (projection[0])[y][i];
		
		Vec3f vtmp2 = cross(vtmp[0], vtmp[1]);
		vtmp2.normalize();
		for (int y = 0; y < 3; ++y)
			ans[y] = vtmp2[y];
		ans[3] = 0.0;
	}
	else {
		Matrix3f A;
		Vec3f b;
		for (int y = 0; y < 3; ++y) 
		{
			for (int x = 0; x < 3; ++x)
				A[x][y] = (projection[0])[x][y];
			b[y] = - (projection[0])[3][y];
		}
		Matrix3f iA;
		iA = A.getInverseMatrix();
		b = iA * b;
		
		for (int y = 0; y < 3; ++y)
			ans[y] = b[y];
		ans[3] = 1.0;
	}
	return ans;
}

// get scale
float PerspectiveCamera::getScale(const Vec4f& coord, const int level) const 
{
	if (maxLevel <= level) 
	{
		cerr << "Level is not within a range: " << level << ' ' << maxLevel << endl;
		exit (1);
	}
	
	// For orthographic case
	if ((projection[0])[0][2] == 0.0 && (projection[0])[1][2] == 0.0 &&
		(projection[0])[2][2] == 0.0) 
	{
		const Vec3f xaxis(projection[0][0][0], projection[0][1][0],
						  projection[0][2][0]);
		const Vec3f yaxis(projection[0][0][1], projection[0][1][1],
						  projection[0][2][1]);
		return (0x0001 << level) / ((xaxis.length() + yaxis.length()) / 2.0);
	}
	else 
	{      
		//const float fz = coord * m_projection[level][2];    
		//return fz * (0x0001 << level) / m_ipscale;
		// ???? new by take into angle difference
		Vec4f ray = coord - center;
		return ray.length() * (0x0001 << level) / ipscale;
	}
}

void PerspectiveCamera::setProjection(std::vector<float>& intrinsics,
						   std::vector<float>& extrinsics,
						   Matrix4f& projection) 
{
	double params[12];
	for (int i = 0; i < 6; ++i) 
	{
		params[i] = intrinsics[i];
		params[6 + i] = extrinsics[i];
	}
	
	for (int y = 0; y < 3; ++y) 
	{
		for (int x = 0; x < 4; ++x ) 
		{
			projection[x][y] = params[4 * y + x];
		}
	}
}

void PerspectiveCamera::setProjectionSub(double params[], Matrix4f& projection, const int level) 
{
	const double rx = params[6] * M_PI / 180.0;
	const double ry = params[7] * M_PI / 180.0;
	const double rz = params[8] * M_PI / 180.0;
	
	const double fovx = params[0] * M_PI / 180.0;
	
	const double f = params[1] / 2.0 / tan(fovx / 2.0);
	Matrix3f K;
	K[0] = Vec3f(f, 0.0, 0.0);
	K[1] = Vec3f(0.0, f, 0.0);
	K[2] = Vec3f(0.0, 0.0, -1.0);
	
	Matrix3f trans;
	trans[0] = Vec3f(1.0, 0.0, 0.0);
	trans[1] = Vec3f(0.0, -1.0, 0.0);
	trans[2] = Vec3f(params[1] / 2.0, params[2] / 2.0, 1.0);
	
	K = trans * K;
	
	Matrix3f Rx;
	Rx[0] = Vec3f(1.0, 0.0, 0.0);
	Rx[1] = Vec3f(0.0f, cos(rx), sin(rx));
	Rx[2] = Vec3f(0.0, -sin(rx), cos(rx));
	
	Matrix3f Ry;
	Ry[0] = Vec3f(cos(ry), 0, -sin(ry));
	Ry[1] = Vec3f(0.0, 1.0, 0.0);
	Ry[2] = Vec3f(sin(ry), 0, cos(ry));
	
	Matrix3f Rz;
	Rz[0] = Vec3f(cos(rz), sin(rz), 0.0);
	Rz[1] = Vec3f(-sin(rz), cos(rz), 0.0);
	Rz[2] = Vec3f(0.0, 0.0, 1.0);
	
	Matrix3f R = Rx.transpose() * Ry.transpose() * Rz.transpose();
	
	Vec3f t(params[3], params[4], params[5]);
	
	Matrix3f left = K * R;
	Vec3f right = - K * (R * t);
	
	for (int y = 0; y < 3; ++y) {
		for (int x = 0; x < 3; ++x)
			projection[x][y] = left[x][y];
		projection[3][y] = right[y];
	}
	
	//halve the first and second rows by 2
	for(int i=0; i<4; i++)
	{
		projection[i][0] /= 2.0;
		projection[i][1] /= 2.0;
	}
	
}

void PerspectiveCamera::projToParam(Matrix4f& mat, double q[6]) {
	double s;
	int i;
	
	q[3] = mat[3][0];
	q[4] = mat[3][1];
	q[5] = mat[3][2];
	q[0] = 0;
	q[1] = 0;
	q[2] = 0;
	if (mat[0][2] == 1.0) 
	{
		q[1] = (double) -M_PI/2.0;
		q[2] = 0;
		q[0]=atan2(-mat[1][0],mat[1][1]);
	}
	else 
	{
		if (mat[0][2] == -1.0) 
		{ 
			q[1] = M_PI/2.0;
			q[2] = 0;
			q[0]=atan2(mat[1][0],mat[1][1]);    
		}
		else 
		{
			q[1] = (double)  asin(-mat[0][2]);
			if (cos(q[1]) > 0.0) 
			{ 
				s = 1.0;
			} 
			else 
			{ 
				s =-1.0;
			}
			q[0] =atan2(mat[1][2]*s, mat[2][2]*s); 
			q[2] =atan2(mat[0][1]*s, mat[0][0]*s); 
		}
	}
	q[0]=q[0]*180/M_PI;//RadInDeg;
	q[1]=q[1]*180/M_PI;//RadInDeg;
	q[2]=q[2]*180/M_PI;//RadInDeg;
	for(i=0;i<3;i++){
		if (fabs(q[i])>180.0){
			q[i]= (q[i]>0) ? q[i]-360.0 : q[i]+360.0;
		}
	}
}

void PerspectiveCamera::paramToProj(const double q[6], Matrix4f& mat) {
	const double a = q[0] * M_PI / 180.0;
	const double b = q[1] * M_PI / 180.0;
	const double g = q[2] * M_PI / 180.0;
	
	const double s1=sin(a);  const double s2=sin(b);  const double s3=sin(g);
	const double c1=cos(a);  const double c2=cos(b);  const double c3=cos(g);
	
	/*   Premiere colonne*/	/*   Seconde colonne	*/
	mat[0][0]=c2*c3; 		mat[1][0]=c3*s2*s1-s3*c1;  
	mat[0][1]=s3*c2; 		mat[1][1]=s3*s2*s1+c3*c1; 
	mat[0][2]=-s2;   		mat[1][2]=c2*s1;
	
	/*   Troisieme colonne*/	/*  Quatrieme colonne	*/
	mat[2][0]=c3*s2*c1+s3*s1; 	mat[3][0]=q[3]; 
	mat[2][1]=s3*s2*c1-c3*s1; 	mat[3][1]=q[4]; 
	mat[2][2]=c2*c1;			mat[3][2]=q[5];
	
	mat[0][3] = mat[1][3] = mat[2][3] = 0.0;
	mat[3][3] = 1.0;
}

float PerspectiveCamera::computeDepthDiff(const Vec4f& lhs, const Vec4f& rhs) const 
{
	// orthographic projection case
	if (projection[0][0][2] == 0.0 && projection[0][1][2] == 0.0 &&
		projection[0][2][2] == 0.0) 
	{
		return - dot(center, (lhs - rhs));
	}
	else 
	{
		return dot(oaxis, (lhs - rhs));
	}
}

float PerspectiveCamera::computeDistance(const Vec4f& point) const 
{
	const float fx = point[0] - center[0];
	const float fy = point[1] - center[1];
	const float fz = point[2] - center[2];
	
	return sqrt(fx * fx + fy * fy + fz * fz);
}

float PerspectiveCamera::computeDepth(const Vec4f& point) const 
{
	// orthographic projection case
	if (projection[0][0][2] == 0.0 && projection[0][1][2] == 0.0 &&
		projection[0][2][2] == 0.0) 
	{
		//cerr << "Because I'm using negative depth to represent flip. this could be a problem" << endl;
		return - dot(center, point);
	}
	else 
	{
		return dot(oaxis, point);
	}
}

void PerspectiveCamera::getPlaneAxes(const Vec4f& coord, const Vec4f& normal,
						  Vec4f& pxaxis, Vec4f& pyaxis, const int level)
{
	const float pscale = getScale(coord, level);
	
	Vec3f normal3(normal[0], normal[1], normal[2]);
	Vec3f yaxis3 = cross(normal3, xaxis);
	yaxis3.normalize();
	Vec3f xaxis3 = cross(yaxis3, normal3);
	pxaxis[0] = xaxis3[0];  pxaxis[1] = xaxis3[1];  pxaxis[2] = xaxis3[2];  pxaxis[3] = 0.0;
	pyaxis[0] = yaxis3[0];  pyaxis[1] = yaxis3[1];  pyaxis[2] = yaxis3[2];  pyaxis[3] = 0.0;
	
	pxaxis *= pscale;
	pyaxis *= pscale;
	const float xdis = (project(coord + pxaxis, level) -
						project(coord, level)).length();
	const float ydis = (project(coord + pyaxis, level) -
						project(coord, level)).length();
	pxaxis *= axesScale / xdis;
	pyaxis *= axesScale / ydis;
	
}

void PerspectiveCamera::setAxesScale(const float axesScale_) 
{
	axesScale = axesScale_;
}  

void PerspectiveCamera::intersect(const Vec4f& coord, const Vec4f& abcd,
					   Vec4f& cross, float& distance) const 
{
	Vec4f ray = coord - center;
	ray.normalize();
	const float A = dot(coord,abcd);
	const float B = dot(ray,abcd);
	
	if (B == 0.0f) {
		distance = 0xffff;
		cross = Vec4f(0.0f, 0.0f, 0.0f, -1.0f);
	}
	else 
	{
		distance = - A / B;
		cross = coord + distance * ray;
	}  
}

Vec4f PerspectiveCamera::intersect(const Vec4f& coord, const Vec4f& abcd) const 
{
	Vec4f ray = center - coord;
	
	const float A = dot(coord,abcd);
	const float B = dot(ray,abcd);
	
	if (B == 0.0f)
		return Vec4f(0.0f, 0.0f, 0.0f, -1.0f);
	else
		return coord - A / B * ray;
}

Vec3f PerspectiveCamera::unprojectPixelwithDepth(const Vec3f& icoord, const int m_level) const 
{
	Matrix3f A;
	Vec3f b(icoord[0], icoord[1], icoord[2]);
	for (int y = 0; y < 3; ++y) {
		for (int x = 0; x < 3; ++x)
			A[y][x] = projection[m_level][y][x];
		b[y] -= projection[m_level][3][y];    
	}
	Matrix3f IA = A.getInverseMatrix();
	Vec3f x = IA * b;
	return x;
}

void PerspectiveCamera::setFundementalMatrix(PerspectiveCamera lhs, PerspectiveCamera rhs, Matrix3f& F, const int level) 
{
	Matrix4f projMatLhs = (lhs.getProjectionMatrix(level)).transpose();
	const Vec4f& p00 = projMatLhs[0];
	const Vec4f& p01 = projMatLhs[1];
	const Vec4f& p02 = projMatLhs[2];
	
	Matrix4f projMatRhs = (rhs.getProjectionMatrix(level)).transpose();
	const Vec4f& p10 = projMatRhs[0];
	const Vec4f& p11 = projMatRhs[1];
	const Vec4f& p12 = projMatRhs[2];
	
	F[0][0] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p01, p02, p11, p12).transpose());
	F[1][0] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p01, p02, p12, p10).transpose());
	F[2][0] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p01, p02, p10, p11).transpose());
	
	F[0][1] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p02, p00, p11, p12).transpose());
	F[1][1] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p02, p00, p12, p10).transpose());
	F[2][1] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p02, p00, p10, p11).transpose());
	
	F[0][2] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p00, p01, p11, p12).transpose());
	F[1][2] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p00, p01, p12, p10).transpose());
	F[2][2] = MathUtils::getDeterminentOfMatrix4f(makeMatrix4f(p00, p01, p10, p11).transpose());
}

float PerspectiveCamera::computeDistanceToEpipolarLine(Matrix3f F, Vec3f p0, Vec3f p1) 
{
	Vec3f line = F * p1;
	const float ftmp = sqrt(line[0] * line[0] + line[1] * line[1]);
	if (ftmp == 0.0)
		return 0.0;
	
	line /= ftmp;
	return fabs(dot(line, p0));
}

void PerspectiveCamera::polarDecompose()
{
    Vec3f a1((projection[0])[0][0], (projection[0])[1][0], (projection[0])[2][0]);
    Vec3f a2((projection[0])[0][1], (projection[0])[1][1], (projection[0])[2][1]);
    Vec3f a3((projection[0])[0][2], (projection[0])[1][2], (projection[0])[2][2]);
	
    double t1 = (projection[0])[3][0], t2 = (projection[0])[3][1], t3 = (projection[0])[3][2];
	
    double rho = 1 / a3.length();
    a1 *= rho; a2 *= rho; a3 *= rho;
    t1 *= rho; t2 *= rho; t3 *= rho;
	
    Vec3f  r3 = a3;
    double u0 = dot(a1,r3);
    double v0 = dot(a2,r3);
	
    Vec3f v = a2 - r3 * v0;
    double beta = v.length();
    Vec3f r2 = v / beta;
	
    double matt = dot(a1,r2);
    v = a1 - r2 * matt - r3 * u0;
    double alpha = v.length();
    Vec3f r1 = v / alpha;
	
    double tz = t3;
    double ty = (t2 - v0 * tz) / beta;
    double tx = (t1 - u0 * tz - ty * matt) / alpha;
	
    intrinsicMat.setIdentity();
	extrinsicMat.setIdentity();
	
    intrinsicMat[0][0] = alpha; intrinsicMat[1][0] = matt;  intrinsicMat[2][0] = u0;
    intrinsicMat[0][1] = 0;     intrinsicMat[1][1] = beta;  intrinsicMat[2][1] = v0;
	
    extrinsicMat[0][0] = r1[0]; extrinsicMat[1][0] = r1[1]; extrinsicMat[2][0] = r1[2]; extrinsicMat[3][0] = tx;
    extrinsicMat[0][1] = r2[0]; extrinsicMat[1][1] = r2[1]; extrinsicMat[2][1] = r2[2]; extrinsicMat[3][1] = ty;
    extrinsicMat[0][2] = r3[0]; extrinsicMat[1][2] = r3[1]; extrinsicMat[2][2] = r3[2]; extrinsicMat[3][2] = tz;

}

void PerspectiveCamera::updateProjectionMatrixWithIntrinsicAndExtrinsic()
{
	setProjectionMatrix(expand3To4(intrinsicMat) * extrinsicMat);
}