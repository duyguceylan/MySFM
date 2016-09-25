/*
 *  Photo.cpp
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "Photo.h"

#include <fstream>

using namespace std;

Photo::Photo(const std::string imgName, const int maxLevel_) 
{
	maxLevel = maxLevel_;
	imgPyramid.resize(maxLevel);
	for(int i=0; i<maxLevel; i++)
		imgPyramid[i] = NULL;
	
	//read image
	imgPyramid[0] = new Img();
	imgPyramid[0]->read(imgName);
	
	//build the image pyramid
	buildImagePyramid();
	
	cam = NULL;
}

Photo::Photo(const std::string imgName, const std::string camName, const int maxLevel_) 
{
	maxLevel = maxLevel_;
	imgPyramid.resize(maxLevel);
	for(int i=0; i<maxLevel; i++)
		imgPyramid[i] = NULL;
	
	//read image
	imgPyramid[0] = new Img();
	imgPyramid[0]->read(imgName);
	//initialize camera
	cam = new PerspectiveCamera();
	cam->init(camName, maxLevel);
	//build the image pyramid
	buildImagePyramid();
}

Photo::~Photo() 
{
	if(cam != NULL)
		delete cam;
	for(int i=0; i<maxLevel; i++)
		delete imgPyramid[i];
	imgPyramid.clear();
}

Color Photo::getColor(const float fx, const float fy, const int level) const 
{
	return imgPyramid[level]->lerp(fx, fy);
};

Color Photo::getColor(const Vec4f& coord, const int level) const 
{
	const Vec3f icoord = cam->project(coord, level);
	return imgPyramid[level]->lerp(icoord[0], icoord[1]);
};

void Photo::buildImagePyramid()
{
	if(imgPyramid[0] == NULL)
	{
		cerr << "Level 0 of image pyramid is NULL, cannot build the pyramid." << endl;
		return;
	}
	
	//apply gaussian filter
	Matrix4f mask;
	mask[0] = Vec4f(1.0, 3.0, 3.0, 1.0);  mask[1] = Vec4f(3.0, 9.0, 9.0, 3.0);
	mask[2] = Vec4f(3.0, 9.0, 9.0, 3.0);  mask[3] = Vec4f(1.0, 3.0, 3.0, 1.0);
	
	float total = 64.0f;
	mask /= total;
	
	//----------------------------------------------------------------------
	// fill the pyramid
	for (int level = 1; level < maxLevel; ++level) 
	{
		int wPrev = imgPyramid[level-1]->width();
		int hPrev = imgPyramid[level-1]->height();
		int w = wPrev / 2;
		int h = hPrev / 2;
		imgPyramid[level] = new Img(w, h);
		for (int y = 0; y < h; ++y) 
		{      
			for (int x = 0; x < w; ++x) 
			{
				Color color(0.0, 0.0, 0.0);
				Color tmpColor;
				float denom = 0.0;
				
				for (int j = -1; j < 3; ++j) 
				{
					int ytmp = 2 * y + j;
					if (ytmp < 0 || hPrev - 1 < ytmp)
						continue;
					
					for (int i = -1; i < 3; ++i) 
					{
						int xtmp = 2 * x + i;
						if (xtmp < 0 ||wPrev - 1 < xtmp)
							continue;
						
						tmpColor = (*imgPyramid[level-1])(xtmp, ytmp);
						tmpColor[0] = mask[i+1][j+1] * tmpColor[0];
						tmpColor[1] = mask[i+1][j+1] * tmpColor[1];
						tmpColor[2] = mask[i+1][j+1] * tmpColor[2];
						
						color += tmpColor;
						
						denom += mask[i+1][j+1];
					}
				}
				
				color /= denom;
				imgPyramid[level]->setColor(x,y,color);
			}
		}
	}
}

float Photo::computeSSD(const std::vector<Color>& tex0, const std::vector<Color>& tex1) 
{
	float ans = 0.0f;
	for (int i = 0; i < (int)tex0.size(); ++i)
		ans += dot(tex0[i] - tex1[i], tex0[i] - tex1[i]);
	
	//max difference: (1.0 * 1.0)*3
	ans /= (int)tex0.size() * 3.0;
	
	return ans;
}

void Photo::normalizeImgPatch(std::vector<Color>& tex) {
	//----------------------------------------------------------------------
	// normalize average
	Vec3f ave(0.0, 0.0, 0.0);
	for (int i = 0; i < (int)tex.size(); ++i)
		ave += tex[i];
	ave /= (int)tex.size();
	
	for (int i = 0; i < (int)tex.size(); ++i)
		tex[i] -= ave;
	
	//----------------------------------------------------------------------  
	// compute variance
	float ave2 = 0.0f;
	for (int i = 0; i < (int)tex.size(); ++i)
		ave2 += dot(tex[i], tex[i]);
	ave2 /= (int)tex.size() * 3;
	ave2 = sqrt(ave2);
	if (ave2 == 0.0f)
		ave2 = 1.0f;
	
	for (int i = 0; i < (int)tex.size(); ++i)
		tex[i] /= ave2;
}

void Photo::grabTexture(const int level, const Vec2f& icoord,
						const Vec2f& xaxis, const Vec2f& yaxis,
						const int size, std::vector<Color>& tex,
						const int normalizef) const
{
	const int margin = size / 2;
	
	// Check boundary condition
	const float maxx = icoord[0] + size * fabs(xaxis[0]) + size * fabs(yaxis[0]);
	const float minx = icoord[0] - size * fabs(xaxis[0]) - size * fabs(yaxis[0]);
	const float maxy = icoord[1] + size * fabs(xaxis[1]) + size * fabs(yaxis[1]);
	const float miny = icoord[1] - size * fabs(xaxis[1]) - size * fabs(yaxis[1]);
	
	tex.clear();
	int w = imgPyramid[level]->width();
	int h = imgPyramid[level]->height();
	
	if (minx < 0 || w - 1 <= maxx ||
		miny < 0 || h - 1 <= maxy)
		return;
	
	for (int y = -margin; y <= margin; ++y) 
	{
		Vec2f v2ftmp = icoord - xaxis*margin + yaxis*y;
		for (int x = -margin; x <= margin; ++x) 
		{
			tex.push_back((*imgPyramid[level])(v2ftmp[0], v2ftmp[1]));
			//imgPyramid[level]->setColor((int)v2ftmp[0], (int)v2ftmp[1], Color(1.0, 0.0, 0.0));
			v2ftmp += xaxis;
		}
	}
	
	if (normalizef)
		normalizeImgPatch(tex);
}

void Photo::grabTexture(const int level, const Vec4f& coord,
						const Vec4f& pxaxis, const Vec4f& pyaxis, const Vec4f& pzaxis,
						const int size, std::vector<Color>& tex, float& weight,
						const int normalizef) const 
{
	const int scale = 0x0001 << level;
	
	const Vec3f icoord3 = cam->project(coord, level);
	const Vec2f icoord(icoord3[0], icoord3[1]);
	
	const Vec3f xaxis3 = cam->project(coord + pxaxis * scale, level) - icoord3;
	const Vec2f xaxis(xaxis3[0], xaxis3[1]);
	
	const Vec3f yaxis3 = cam->project(coord + pyaxis * scale, level) - icoord3;
	const Vec2f yaxis(yaxis3[0], yaxis3[1]);
	
	grabTexture(level, icoord, xaxis, yaxis, size, tex, normalizef);
	
	Vec4f ray = cam->getCenter() - coord;
	ray.normalize();
	weight = max(0.0f, dot(pzaxis,ray));
}