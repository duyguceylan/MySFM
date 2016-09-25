//
//  OpenGLViewPort.h
//  BundlerTrackViewer
//
//  Created by Duygu Ceylan on 11/8/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#ifdef __APPLE_CC__
#include "OpenGL/gl.h"
#include "OpenGL/glu.h"
#endif

#include <vector>
#include "Image.h"
#include "Common.h"

class BundlerOutputParser;

enum ImageType
{
	HORIZONTAL = 0,
	VERTICAL
};

@interface OpenGLViewPort : NSOpenGLView {
	
	unsigned int openGLScreenWidth,openGLScreenHeight;
	
	GLint viewportDimensions[4];
	float textureIntensity;
	float backgroundColor[4];
	int numberOfScreenshotsTaken;
	
	BundlerOutputParser *bundleParser;
	
	std::vector<std::string> listOfInputImagePaths;
	std::vector<float * > allInputImages;
	bool allImagesLoaded;
	
	bool grids3DComputed;
	bool drawLine;
	
	bool showRectifiedImage;
	float *rectifiedImage;
	int rectImageWidth, rectImageHeight;
	float rectImageScale;
	GLuint rectifiedTexture;
	
	int currentImageIndex;
	int currentImageWidth, currentImageHeight;
	float * currentImageBitmap;
	GLuint currentTexture;
	Vec2f currentFeature;
	vector<Vec2f> currentFeatures;
	
	int linkedImageIndex;
	int linkedImageWidth, linkedImageHeight;
	float * linkedImageBitmap;
	GLuint linkedTexture;
	Vec2f linkedFeature;
	vector<Vec2f> linkedFeatures;
	bool isValidCorr;
	vector<bool> validCorr;
	
	vector<Vec2f> epipolarLine;
	
	vector<vector<Vec2f> > currentRepetitions;
	vector<Color> currentRepetitionColors;
	vector<vector<Vec2f> > currentRefinedRepetitions;
	
	keypt_t *currentFeatureInfo;
	int currentFeatureNo;
	
	float *pointProjections;
	Color *pointProjectionColors;
	int noProjectedPoints;
	
	vector<vector<Vec2f> > currentGridStrokes;
	
	ImageType type;
	float scaleOfImg;
	
	float imgZoomScale;
	Vec2f imgTranslate;
	float maxImgTranslate;
	float imgTranslateDelta;
	float txmin, tymin, txmax, tymax; //texture coordinates
	float sxmin, symin, sxmax, symax; //screen coordinates
	Vec2f currentCoord;
	NSPoint startLocation;
	vector<Vec2f> lineScreenCoordinates;
}

- (void) setScaleOfImg:(float)_scale;
- (void) setFeaturePointInfo:(keypt_t*)featurePointInfo 
			   withFeatureNo:(int)noOfFeatures;
- (void) setProjectedPointInfo:(float*)projectedPoints
					 withColor: (Color*)projectedPointColors
				  withPointNo:(int)noOfPoints;

- (void) setRepetitions: (vector<vector<Vec2f> >)_repetitions
			 withColors:(vector<Color> )_colors;

-(void) setRefinedRepetitions:(vector<vector<Vec2f> >)_refinedRepetitions;

- (void) setGridStrokes: (vector<vector<Vec2f> >) _gridStrokes;

- (void) takeScreenshot;
- (void) takeScreenshot:(NSString *) _fileName;
- (void) renderOffScreenWithWidth:(unsigned int) width andHeight: (unsigned int) height andData: (float *) data;
- (void) setTextureIntensity:(float) v;

- (void)display;
- (void)makeCurrent;
- (void)reshape;
- (void)drawRect:(NSRect)r;

-(void) attachBundlerParser: (BundlerOutputParser *) _bundleParser;
-(void) addNewImagePath: (std::string) filename;
-(void) updateRectifiedImage: (Img*) _rectifiedImage;
-(void) updateActiveImage: (int) _currentImageIndex;
-(void) updateActiveImages: (int) _currentImageIndex
				  linkedTo: (int) _linkedImageIndex;
-(void) updateActiveRectImages: (Img *) _currentImage 
					  linkedTo: (Img *) _linkedImage;
-(void)updateActiveCorrespondences:(Vec2f) _currentFeature
						  linkedTo:(Vec2f) _linkedFeature
						   isValid:(bool) _valid;
-(void)updateActiveCorrespondences:(Vec2f) _currentFeature
						  linkedTo:(Vec2f) _linkedFeature
				 withEpipolarStart:(Vec2f) _start
				   withEpipolarEnd:(Vec2f) _end;
-(void)updateCorrespondences:(vector<Vec2f>) _currentFeatures
					linkedTo:(vector<Vec2f>) _linkedFeatures
					validity:(vector<bool>) _validCorr;

-(void) loadImage: (int) imageIndex
		withWidth:(int &) w
	   withHeight:(int &) h
		 withData:(float**) data
		withScale:(float) scale;

-(void) loadAllImages: (float) scale;

-(void) translateImageLeft;
-(void) translateImageRight;
-(void) translateImageUp;
-(void) translateImageDown;
-(void) scaleImage;
-(void) setDrawLine:(bool) _drawLine;
-(void) setGrids3DComputed:(bool) _grids3DComputed;

-(void) mouseDown:(NSEvent *)theEvent;
-(void) mouseDragged:(NSEvent *)theEvent;
-(void) mouseUp:(NSEvent *)theEvent;
-(void) getTemplateCoord:(Vec2f&) topLeft
				   width:(int&)w
				  height:(int&)h;
-(void) getStrokes:(vector<Vec2f>&) endpoints;
-(void) getCurrentCoord:(Vec2f &)c;

@end
