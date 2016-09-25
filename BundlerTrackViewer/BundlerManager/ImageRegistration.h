/*
 *  ImageRegistration.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 3/7/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _IMAGE_REGISTRATION_H
#define _IMAGE_REGISTRATION_H

class ImageTransformation
{
public:
	int getImageIndex1() {return imageIndex1; }
	int getImageIndex2() {return imageIndex2; }
	int getShiftInColumns() {return shiftCol; }
	int getShiftInRows() {return shiftRow; }
	float getRectificationScale() {return scale; }
	Vec2f getRectificationTranslation() {return translation; }
	Matrix3f getFundementalMatrix() {return fundMatrix; }
	
	void setImageIndex1(int imageIndex1_) {imageIndex1 = imageIndex1_; }
	void setImageIndex2(int imageIndex2_) {imageIndex2 = imageIndex2_; }
	void setShiftInColumns(int shiftCol_ ) { shiftCol = shiftCol_;}
	
private:
	int imageIndex1;
	int imageIndex2;
	
	//for rectified grid-based registration
	int shiftCol;
	int shiftRow;
	float scale;
	Vec2f translation;
	
	//for normal registration
	Matrix3f fundMatrix;
	
	float score;
	int group;
	bool used;
}

#endif