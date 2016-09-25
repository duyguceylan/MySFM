/*
 *  Photo.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _PHOTO_H
#define _PHOTO_H

#include "Image.h"
#include "PerspectiveCamera.h"

class Photo
{
public:
	Photo(const std::string imgName, const int maxLevel_ = 1);
	Photo(const std::string imgName, const std::string camName, const int maxLevel = 1);
	~Photo();
	
	PerspectiveCamera* getCam() {return cam; }
	Img* getImg(int level) {return imgPyramid[level]; }
	float* getColors(int level) {return imgPyramid[level]->getColors(); }
	
	void buildImagePyramid();
	
	// grabTexture centered at the given pixel location
	void grabTexture(const int level, const Vec2f& icoord,
					 const Vec2f& xaxis, const Vec2f& yaxis, const int size,
					 std::vector<Color>& tex, const int normalizef = 1) const;
	
	// grabTexture centered at the projection of the given 3D point
	void grabTexture(const int level, const Vec4f& coord,
					 const Vec4f& pxaxis, const Vec4f& pyaxis, const Vec4f& pzaxis,
					 const int size, std::vector<Color>& tex, float& weight,
					 const int normalizef = 1) const;
	
	Color getColor(const float fx, const float fy, const int level) const;  
	Color getColor(const Vec4f& coord, const int level) const;
	Color getColor(const int x, const int y, const int level) {return (*imgPyramid[level])(x,y);}
	
	void setColor(int x, int y, Color c, int level) {imgPyramid[level]->setColor(x,y,c); }
	
	void saveImage(int level, const string& filename) {imgPyramid[level]->write(filename); }
	
	static void normalizeImgPatch(std::vector<Color>& tex);
	
	static float computeSSD(const std::vector<Color>& tex0, const std::vector<Color>& tex1);
	
	int getWidth(int level) {return imgPyramid[level]->width(); }
	int getHeight(int level) {return imgPyramid[level]->height(); }
	
protected:
	PerspectiveCamera* cam;
	std::vector<Img*> imgPyramid;
	int maxLevel;
};

#endif // PHOTO_H
