/*
 *  Line.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 6/29/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _LINE_H
#define _LINE_H

#include "../MathUtils/MyMatrix.h"
#include "Image.h"

#include <vector>
#include "opencv2/core/core_c.h"

using namespace std;

typedef enum LINERELATION
{
	PARALLEL = 0,
	SAME,
	INTERSECTING,
	NONCOPLANAR,
	NONINTERSECTING
};

typedef struct MatchingLine
{
	int imgIndex;
	int lineIndex;
	float matchingScore;
	Vec2f startPos;
	Vec2f endPos;
	
	MatchingLine(int imgIndex_, int lineIndex_, float matchingScore_) {imgIndex = imgIndex_; lineIndex = lineIndex_; matchingScore = matchingScore_;}
	
	bool operator<(const MatchingLine& rhs) const
	{
		return matchingScore > rhs.matchingScore;
	}
};

typedef struct ImgProjection
{
	int imgIndex;
	Vec2f startPos;
	Vec2f endPos;
	
	ImgProjection(int imgIndex_, Vec2f startPos_, Vec2f endPos_) {imgIndex = imgIndex_; startPos = startPos_; endPos = endPos_; }
};

class Line
{
protected:	
	Vec3f startPoint3D, endPoint3D;
	Vec2f startPoint2D, endPoint2D;
	Vec2f clippedStartPoint, clippedEndPoint;
	Vec2f normal2D;
	float distance2D;
	
	bool matched3DLine;
	bool isContour;
	
	vector<MatchingLine> matches;
	vector<ImgProjection> imageProjections;
	
public:
	Line();
	Line(Vec3f normal, float distance);
	Line(Vec3f startPoint_, Vec3f endpoint_);
	Line(Vec2f projStarPoint_, Vec2f projEndPoint_);
	
	Vec3f getStartPoint3D() {return startPoint3D; }
	Vec3f getEndPoint3D() {return endPoint3D; }
	Vec2f getStartPoint2D() {return startPoint2D; }
	Vec2f getEndPoint2D() {return endPoint2D; }
	Vec2f getClippedStartPoint() {return clippedStartPoint; }
	Vec2f getClippedEndPoint() {return clippedEndPoint; }
	bool isMatched3DLine() {return matched3DLine; }
	bool isContourLine() {return isContour; }
	int getNoMatchingLines() {return matches.size(); }
	MatchingLine getMatchingLine(int index) {return matches[index]; }
	int getNoImgProjections() {return imageProjections.size(); }
	ImgProjection getImgProjection(int index) {return imageProjections[index]; }
	int doesImageHaveProjection(int index);
	
	void setStartPoint3D(Vec3f startPoint3D_) {startPoint3D = startPoint3D_; }
	void setEndPoint3D(Vec3f endPoint3D_) {endPoint3D = endPoint3D_; }
	void setStartPoint2D(Vec2f startPoint2D_) {startPoint2D = startPoint2D_; }
	void setEndPoint2D(Vec2f endPoint2D_) {endPoint2D = endPoint2D_; }
	void setClippedStartPoint(Vec2f clippedStartPoint_) {clippedStartPoint = clippedStartPoint_; }
	void setClippedEndPoint(Vec2f clippedEndPoint_) {clippedEndPoint = clippedEndPoint_;}
	void setMatched3DLine(bool matched3DLine_) {matched3DLine = matched3DLine_; }
	void setContourLine(bool isContour_) {isContour = isContour_; }
	void addNewMatchingLine(MatchingLine l) {matches.push_back(l); }
	void clearMatchingLines() {matches.clear(); }
	void addNewImgProjection(ImgProjection p) {imageProjections.push_back(p);}
	
	Vec2f generateRandomPointOnLine2D();
	float computeLineLength2D();
	Vec2f computeLineNormal2D();
	float computeLineDistance2D();
	void findBoundingRectangle(vector<Vec2f> &boundingRect);
	
	float findPerpendicularDistanceTo2DLine(Vec2f pixelPos);
	void projectToImage(Matrix4f projMatrix);
	Vec3f compute2DLineNormalFormParams();
	
	bool isVisibleByCamera(Matrix4f projMatrix, Vec2f &projStart, Vec2f &projEnd, int width, int height);
	bool isInsideBoundingRectangle(vector<Vec2f> &boundingRect, Vec2f pt);
	
	LINERELATION findRelationWith3DLine(Line *line, Vec3f &intersectionPt);
	LINERELATION findRelationWith2DLine(Line *line, Vec2f &intersectionPt, bool extend = false);
	
	void fillLineMask(int width, int height, Img *tmp);
	
};

#endif