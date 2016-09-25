/*
 *  Camera.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _PERSPECTIVE_CAMERA_H
#define _PERSPECTIVE_CAMERA_H

#include <vector>
#include <string>
#include <climits>
#include "../MathUtils/MyMatrix.h"

class PerspectiveCamera 
{
public:
	PerspectiveCamera(void);
	virtual ~PerspectiveCamera();
	
	void setProjectionMatrix(Matrix4f mat);
	
	// Update projection matrices from intrinsics and extrinsics
	void updateProjection(void);
	// Update all the camera related parameters
	void updateCamera(void);
	
	void init(const std::string cname, const int maxLevel);
	void write(const std::string file);
	
	inline Vec3f project(Vec4f coord, int level);
	inline Vec3f mult(Vec4f coord, int level);
	
	static void setProjection(std::vector<float>& intrinsics,
							  std::vector<float>& extrinsics,
							  Matrix4f& projection);
	
	float getScale(const Vec4f& coord, const int level) const;
	void getPlaneAxes(const Vec4f& coord, const Vec4f& normal,
					  Vec4f& pxaxis, Vec4f& pyaxis, const int level = 0);
	
	void setAxesScale(const float axesScale);
	
	static void projToParam(Matrix4f& mat, double q[6]); //get rotation and translation params
	static void paramToProj(const double q[6], Matrix4f& mat);  //construct from rotation and translation params
	static void setProjectionSub(double params[], Matrix4f& projection,
								 const int level);
	
	float computeDistance(const Vec4f& point) const;
	float computeDepth(const Vec4f& point) const;  
	float computeDepthDiff(const Vec4f& rhs, const Vec4f& lhs) const;
	
	// Compute where the viewing ray passing through coord intersects
	// with the plane abcd.
	Vec4f intersect(const Vec4f& coord, const Vec4f& abcd) const;
	void intersect(const Vec4f& coord, const Vec4f& abcd,
				   Vec4f& cross, float& distance) const;
	
	// Compute the 3D ray starting at the camera position and passing through the given pixel
	Vec3f unprojectPixelwithDepth(const Vec3f& icoord, const int m_level) const;
	
	std::string getCamParamFilename() {return camParamFilename; }
	Vec4f getCenter() {return center; }
	Vec4f getOAxis() {return oaxis; }
	
	void setCamParamFilename(std::string name_) {camParamFilename = name_; }
	void setCenter(Vec4f c_) {center = c_; }
	void setOAxis(Vec4f oa_) {oaxis = oa_; }
	
	Matrix4f getProjectionMatrix(int level) {return projection[level]; }
	float getElementOfProjectionMatrix(int level, int row, int column) {return projection[level][column][row]; }
	
	Matrix4f getExtrinsicMatrix() {return extrinsicMat; }
	Matrix3f getIntrinsicMatrix() {return intrinsicMat; }
	
	void setExtrinsicMatrix(Matrix4f _extrinsicMat) {extrinsicMat = _extrinsicMat; }
	static void setFundementalMatrix(PerspectiveCamera lhs, PerspectiveCamera rhs, Matrix3f& F, const int level = 0);
	static float computeDistanceToEpipolarLine(Matrix3f F, Vec3f p0, Vec3f p1);
	
	void polarDecompose();
	void updateExtrinsicMatrix(Matrix4f trans) {extrinsicMat = trans*extrinsicMat;}
	void scaleIntrinsicMatrix(float scale) {intrinsicMat[0][0] = intrinsicMat[0][0]*scale; intrinsicMat[1][1] = intrinsicMat[1][1]*scale;}
	void updateProjectionMatrixWithIntrinsicAndExtrinsic();
	
protected:
	// txt file name
	std::string camParamFilename;  
	// Optical center
	Vec4f center;
	// Optical axis
	Vec4f oaxis;
	
	float ipscale;
	// 3x4 projection matrix
	std::vector<Matrix4f > projection;
	Vec3f xaxis;
	Vec3f yaxis;
	Vec3f zaxis;
	
	// intrinsic and extrinsic camera parameters. Compact form.
	std::vector<float> intrinsics;
	std::vector<float> extrinsics;
	Matrix3f intrinsicMat;
	Matrix4f extrinsicMat;
	
	int maxLevel;
	
	float axesScale;
	
	Vec4f getOpticalCenter(void) const;
};

inline Vec3f PerspectiveCamera::project(Vec4f coord, int level) 
{
	Vec3f vtmp;    
	
	vtmp = shrink4To3(projection[level] * coord);
	if (vtmp[2] <= 0.0) 
	{
		vtmp[0] = -0xffff;
		vtmp[1] = -0xffff;
		vtmp[2] = -1.0f;
		return vtmp;
	}
	else
		vtmp /= vtmp[2];
	
	vtmp[0] = std::max((float)(INT_MIN + 3.0f),
					   std::min((float)(INT_MAX - 3.0f), vtmp[0]));
	vtmp[1] = std::max((float)(INT_MIN + 3.0f),
					   std::min((float)(INT_MAX - 3.0f), vtmp[1]));
	return vtmp;
};

inline Vec3f PerspectiveCamera::mult(Vec4f coord, int level)
{
	Vec3f vtmp;    
	
	vtmp = shrink4To3(projection[level] * coord);
	
	return vtmp;
};

#endif // PERSPECTIVE_CAMERA_H
