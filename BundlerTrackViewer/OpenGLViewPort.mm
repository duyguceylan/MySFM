//
//  OpenGLViewPort.mm
//  BundlerTrackViewer
//
//  Created by Duygu Ceylan on 11/8/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "OpenGLViewPort.h"
#include "Parser/BundlerOutputParser.h"
#include "BundlerManager/Image.h"

@implementation OpenGLViewPort

- (id)initWithFrame:(NSRect)frameRect
		pixelFormat:(NSOpenGLPixelFormat *) pixFmt
{
	self = [super initWithFrame:frameRect pixelFormat:pixFmt];
	[self setWantsBestResolutionOpenGLSurface:YES];
    
	currentImageIndex = -1;
	linkedImageIndex = -1;
	numberOfScreenshotsTaken = 0;
	type = VERTICAL;
	currentImageBitmap = NULL;
	linkedImageBitmap = NULL;
	rectifiedImage = NULL;
	allImagesLoaded = false;
	currentFeature = Vec2f(0.0, 0.0);
	linkedFeature = Vec2f(0.0, 0.0);
	currentFeatureInfo = NULL;
	
	epipolarLine.clear();
	
	imgZoomScale = 1.0;
	imgTranslate = Vec2f(0.0, 0.0);
	
	return self;
}

- (id)initWithCoder:(NSCoder *)coder
{	
	
	[super initWithCoder: coder];
	
	showRectifiedImage = false;
	rectifiedImage = NULL;
	drawLine = false;
	grids3DComputed = false;
	currentImageIndex = -1;
	linkedImageIndex = -1;
	numberOfScreenshotsTaken = 0;
	type = VERTICAL;
	currentImageBitmap = NULL;
	linkedImageBitmap = NULL;
	allImagesLoaded = false;
	currentFeature = Vec2f(0.0, 0.0);
	linkedFeature = Vec2f(0.0, 0.0);
	currentFeatureInfo = NULL;
	currentRepetitions.clear();
	pointProjections = NULL;
	pointProjectionColors = NULL;
	noProjectedPoints = 0;
	isValidCorr = true;
	currentGridStrokes.clear();
	
	imgZoomScale = 1.0;
	imgTranslate = Vec2f(0.0, 0.0);
	maxImgTranslate = (1.0 - imgZoomScale)/2.0;
	imgTranslateDelta = maxImgTranslate / 5.0;
	
	NSOpenGLPixelFormatAttribute attributes [] = {	NSOpenGLPFAWindow,
													NSOpenGLPFADepthSize, 
													(NSOpenGLPixelFormatAttribute)32,
													NSOpenGLPFASampleBuffers, 
													(NSOpenGLPixelFormatAttribute)1,
													NSOpenGLPFASamples, 
													(NSOpenGLPixelFormatAttribute)24,
													NSOpenGLPFANoRecovery,
													(NSOpenGLPixelFormatAttribute)nil};
	
	NSOpenGLPixelFormat* pixFmt = [[NSOpenGLPixelFormat alloc] initWithAttributes: attributes] ;
	
	[self setPixelFormat:pixFmt];
	[self setWantsBestResolutionOpenGLSurface:YES];
    
	glGenTextures( 1, &currentTexture );
	glGenTextures(1, &linkedTexture);
	glGenTextures(1, &rectifiedTexture);
	
	numberOfScreenshotsTaken = 0;
	type = VERTICAL;
	allImagesLoaded = false;
	
	backgroundColor[0] = 0.2;
	backgroundColor[1] = 0.2;
	backgroundColor[2] = 0.2;
	backgroundColor[3] = 1.0;
	
	textureIntensity = 1.0;
	
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	return self;
}

-(void)display
{
	GLint VBL = 1;								
	[[self openGLContext] setValues:&VBL forParameter:NSOpenGLCPSwapInterval];
	
	[self makeCurrent];
	
	glClearColor(backgroundColor[0],backgroundColor[1],backgroundColor[2],backgroundColor[3]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glViewport(0,0,openGLScreenWidth,openGLScreenHeight);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
    glLoadIdentity();
    glOrtho(0, openGLScreenWidth, openGLScreenHeight, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
    glLoadIdentity();	
	
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	
	if(showRectifiedImage)
	{
		glEnable(GL_COLOR_MATERIAL);
		glEnable( GL_TEXTURE_2D );
		
		// set up texture
		glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
		// when texture area is small, bilinear filter the closest mipmap
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
						GL_LINEAR_MIPMAP_NEAREST );
		// when texture area is large, bilinear filter the original
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
		
		// the texture wraps over at the edges (repeat)
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
		
		gluBuild2DMipmaps( GL_TEXTURE_2D, 4, rectImageWidth, rectImageHeight,
						  GL_RGBA, GL_FLOAT, rectifiedImage);
		
		glBindTexture( GL_TEXTURE_2D, rectifiedTexture );
		
		glColor4f(textureIntensity,textureIntensity,textureIntensity, 1.0f);
		
		int width = floor(rectImageWidth/rectImageScale);
		int height = floor(rectImageHeight/rectImageScale);
		
		txmin = (1.0-imgZoomScale)/2.0;
		txmax = txmin + imgZoomScale;
		tymin = (1.0-imgZoomScale)/2.0;
		tymax = tymin + imgZoomScale;
		txmin += imgTranslate[0];
		txmax += imgTranslate[0];
		tymin += imgTranslate[1];
		tymax += imgTranslate[1];
		
		glBegin(GL_QUADS);
		{
			glTexCoord2f(txmin,tymin);
			glVertex2f((openGLScreenWidth-width)/2.0,(openGLScreenHeight-height)/2.0);
			glTexCoord2f(txmin,tymax);
			glVertex2f((openGLScreenWidth-width)/2.0,(openGLScreenHeight-height)/2.0+height);
			glTexCoord2f(txmax,tymax);
			glVertex2f((openGLScreenWidth-width)/2.0+width,(openGLScreenHeight-height)/2.0+height);
			glTexCoord2f(txmax,tymin);
			glVertex2f((openGLScreenWidth-width)/2.0+width,(openGLScreenHeight-height)/2.0);
		}
		glEnd();
		glDisable( GL_TEXTURE_2D );
		glDisable(GL_COLOR_MATERIAL);

		glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
		glLineWidth(2.0);
		glBegin(GL_LINES);
		for(int i=0; i<lineScreenCoordinates.size(); i+=2)
		{
			glVertex2f(lineScreenCoordinates[i][0], lineScreenCoordinates[i][1]);
			glVertex2f(lineScreenCoordinates[i+1][0], lineScreenCoordinates[i+1][1]);
		}
		glEnd(); 
			
		//repetitions
		glEnable (GL_BLEND);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		for(int i=0; i<currentRepetitions.size(); i++)
		{
			glColor4f(currentRepetitionColors[i][0], currentRepetitionColors[i][1], currentRepetitionColors[i][2], 0.5f);
			Vec2f start((openGLScreenWidth - rectImageWidth/rectImageScale)/2.0, (openGLScreenHeight-rectImageHeight/rectImageScale)/2.0);
			Vec2f end(start[0]+rectImageWidth/rectImageScale, start[1]+rectImageHeight/rectImageScale);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for(int j=0; j<currentRepetitions[i].size(); j++)
			{
				
				Vec2f pos = currentRepetitions[i][j];
				float tX = pos[0] / rectImageWidth; float tY = pos[1] / rectImageHeight;
				float ratioX = (tX - txmin) / (txmax-txmin);
				float ratioY = (tY - tymin) / (tymax-tymin);
				pos[0] = ratioX * (end[0]-start[0]) + start[0];
				pos[1] = ratioY * (end[1]-start[1]) + start[1];
				
				glVertex2f(pos[0],pos[1]);
			}
			glEnd();
		}
		
		glDisable(GL_BLEND);
		
		glFinish();
		return;
	}
	
	if(currentImageBitmap != NULL)
	{
		glEnable(GL_COLOR_MATERIAL);
		glEnable( GL_TEXTURE_2D );
		
		// set up texture
		glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
		// when texture area is small, bilinear filter the closest mipmap
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
						GL_LINEAR_MIPMAP_NEAREST );
		// when texture area is large, bilinear filter the original
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
		
		// the texture wraps over at the edges (repeat)
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
		
		gluBuild2DMipmaps( GL_TEXTURE_2D, 4, currentImageWidth, currentImageHeight,
							  GL_RGBA, GL_FLOAT, currentImageBitmap);
		
		glBindTexture( GL_TEXTURE_2D, currentTexture );
		
		glColor4f(textureIntensity,textureIntensity,textureIntensity, 1.0f);
		
		if(type == VERTICAL)
		{
			if(linkedImageBitmap != NULL)
			{
				glBegin(GL_QUADS);
				{
					glTexCoord2f(0.0,0.0);
					glVertex2f(0.0,openGLScreenHeight/6.0);
					glTexCoord2f(0.0,1.0);
					glVertex2f(0.0,5.0*openGLScreenHeight/6.0);
					glTexCoord2f(1.0,1.0);
					glVertex2f(openGLScreenWidth/2.0,5.0*openGLScreenHeight/6.0);
					glTexCoord2f(1.0,0.0);
					glVertex2f(openGLScreenWidth/2.0,openGLScreenHeight/6.0);
				}
				glEnd();
				glDisable( GL_TEXTURE_2D );
				glDisable(GL_COLOR_MATERIAL);
			}
			else
			{
				txmin = (1.0-imgZoomScale)/2.0;
				txmax = txmin + imgZoomScale;
				tymin = (1.0-imgZoomScale)/2.0;
				tymax = tymin + imgZoomScale;
				txmin += imgTranslate[0];
				txmax += imgTranslate[0];
				tymin += imgTranslate[1];
				tymax += imgTranslate[1];
				
				sxmin = openGLScreenWidth*0.25;
				symin = openGLScreenHeight/6.0;
				sxmax = openGLScreenWidth*0.75;
				symax = 5.0*openGLScreenHeight/6.0;
				
				glBegin(GL_QUADS);
				{
					glTexCoord2f(txmin,tymin);
					glVertex2f(sxmin,symin);
					glTexCoord2f(txmin,tymax);
					glVertex2f(sxmin,symax);
					glTexCoord2f(txmax,tymax);
					glVertex2f(sxmax,symax);
					glTexCoord2f(txmax,tymin);
					glVertex2f(sxmax,symin);
				}
				glEnd();
				glDisable( GL_TEXTURE_2D );
				glDisable(GL_COLOR_MATERIAL);
				
				glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
				glLineWidth(2.0);
				glBegin(GL_LINES);
				for(int i=0; i<lineScreenCoordinates.size(); i+=2)
				{
					glVertex2f(lineScreenCoordinates[i][0], lineScreenCoordinates[i][1]);
					glVertex2f(lineScreenCoordinates[i+1][0], lineScreenCoordinates[i+1][1]);
				}
				glEnd(); 
				
				glEnable(GL_BLEND);
				glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				glColor4f(0.0f, 0.0f, 1.0f, 0.5f);
				glLineWidth(2.0);
				glPointSize(1.0);
				for(int i=0; i<currentGridStrokes.size(); i++)
				{
					glBegin(GL_LINES);
					for(int j=0; j<currentGridStrokes[i].size(); j++)
					{
						int next = j+1;
						if(next == currentGridStrokes[i].size())
							next = 0;
						
						Vec2f c = currentGridStrokes[i][j];
						c /= scaleOfImg;
						c[0] = sxmin + ((sxmax-sxmin)*(c[0]-currentImageWidth*txmin)/(currentImageWidth*(txmax-txmin)));
						c[1] = symin + ((symax-symin)*(c[1]-currentImageHeight*tymin)/(currentImageHeight*(tymax-tymin)));
						if(c[0]>=sxmin && c[0]<=sxmax && c[1]>=symin && c[1]<=symax)
							glVertex2f(c[0],c[1]);
						
						c = currentGridStrokes[i][next];
						c /= scaleOfImg;
						c[0] = sxmin + ((sxmax-sxmin)*(c[0]-currentImageWidth*txmin)/(currentImageWidth*(txmax-txmin)));
						c[1] = symin + ((symax-symin)*(c[1]-currentImageHeight*tymin)/(currentImageHeight*(tymax-tymin)));
						if(c[0]>=sxmin && c[0]<=sxmax && c[1]>=symin && c[1]<=symax)
							glVertex2f(c[0],c[1]);
					}
				}
				glEnd();
				glDisable(GL_BLEND);
				
			}
		}
		else if(type == HORIZONTAL)
		{
			if(linkedImageBitmap != NULL)
			{
				glBegin(GL_QUADS);
				{
					glTexCoord2f(0.0,0.0);
					glVertex2f(openGLScreenWidth/6.0,0.0);
					glTexCoord2f(0.0,1.0);
					glVertex2f(openGLScreenWidth/6.0,openGLScreenHeight/2.0);
					glTexCoord2f(1.0,1.0);
					glVertex2f(5.0*openGLScreenWidth/6.0,openGLScreenHeight/2.0);
					glTexCoord2f(1.0,0.0);
					glVertex2f(5.0*openGLScreenWidth/6.0,0.0);
				}
				glEnd();
				glDisable( GL_TEXTURE_2D );
				glDisable(GL_COLOR_MATERIAL);
			}
			else 
			{
				txmin = (1.0-imgZoomScale)/2.0;
				txmax = txmin + imgZoomScale;
				tymin = (1.0-imgZoomScale)/2.0;
				tymax = tymin + imgZoomScale;
				txmin += imgTranslate[0];
				txmax += imgTranslate[0];
				tymin += imgTranslate[1];
				tymax += imgTranslate[1];
				
				sxmin = openGLScreenWidth/6.0;
				symin = openGLScreenHeight*0.25;
				sxmax = 5.0*openGLScreenWidth/6.0;
				symax = openGLScreenHeight*0.75;
				
				glBegin(GL_QUADS);
				{
					glTexCoord2f(txmin,tymin);
					glVertex2f(sxmin,symin);
					glTexCoord2f(txmin,tymax);
					glVertex2f(sxmin,symax);
					glTexCoord2f(txmax,tymax);
					glVertex2f(sxmax,symax);
					glTexCoord2f(txmax,tymin);
					glVertex2f(sxmax,symin);
				}
				glEnd();
				glDisable( GL_TEXTURE_2D );
				glDisable(GL_COLOR_MATERIAL);
				
				glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
				glLineWidth(2.0);
				glBegin(GL_LINES);
				for(int i=0; i<lineScreenCoordinates.size(); i+=2)
				{
					glVertex2f(lineScreenCoordinates[i][0], lineScreenCoordinates[i][1]);
					glVertex2f(lineScreenCoordinates[i+1][0], lineScreenCoordinates[i+1][1]);
				}
				glEnd(); 
				
				glEnable(GL_BLEND);
				glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				glColor4f(0.0f, 0.0f, 1.0f, 0.5f);
				glLineWidth(2.0);
				glPointSize(1.0);
				
				for(int i=0; i<currentGridStrokes.size(); i++)
				{
					glBegin(GL_LINES);
					for(int j=0; j<currentGridStrokes[i].size(); j++)
					{
						int next = j+1;
						if(next == currentGridStrokes[i].size())
							next = 0;
						
						Vec2f c = currentGridStrokes[i][j];
						c /= scaleOfImg;
						c[0] = sxmin + ((sxmax-sxmin)*(c[0]-currentImageWidth*txmin)/(currentImageWidth*(txmax-txmin)));
						c[1] = symin + ((symax-symin)*(c[1]-currentImageHeight*tymin)/(currentImageHeight*(tymax-tymin)));
						if(c[0]>=sxmin && c[0]<=sxmax && c[1]>=symin && c[1]<=symax)
							glVertex2f(c[0],c[1]);
						
						c = currentGridStrokes[i][next];
						c /= scaleOfImg;
						c[0] = sxmin + ((sxmax-sxmin)*(c[0]-currentImageWidth*txmin)/(currentImageWidth*(txmax-txmin)));
						c[1] = symin + ((symax-symin)*(c[1]-currentImageHeight*tymin)/(currentImageHeight*(tymax-tymin)));
						if(c[0]>=sxmin && c[0]<=sxmax && c[1]>=symin && c[1]<=symax)
							glVertex2f(c[0],c[1]);
					}
					glEnd();
				}
				glDisable(GL_BLEND);
			}
		}
	}
	
	if(linkedImageBitmap != NULL)
	{
		glEnable(GL_COLOR_MATERIAL);
		glEnable( GL_TEXTURE_2D );
		
		// set up texture
		glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
		// when texture area is small, bilinear filter the closest mipmap
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
						GL_LINEAR_MIPMAP_NEAREST );
		// when texture area is large, bilinear filter the original
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
		
		// the texture wraps over at the edges (repeat)
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
		
		gluBuild2DMipmaps( GL_TEXTURE_2D, 4, linkedImageWidth, linkedImageHeight,
						  GL_RGBA, GL_FLOAT, linkedImageBitmap);
		
		glBindTexture( GL_TEXTURE_2D, linkedTexture );
		
		glColor4f(textureIntensity,textureIntensity,textureIntensity, 1.0f);
		
		if(type == VERTICAL)
		{
			glBegin(GL_QUADS);
			{
				glTexCoord2f(0.0,0.0);
				glVertex2f(openGLScreenWidth/2.0,openGLScreenHeight/6.0);
				glTexCoord2f(0.0,1.0);
				glVertex2f(openGLScreenWidth/2.0,5.0*openGLScreenHeight/6.0);
				glTexCoord2f(1.0,1.0);
				glVertex2f(openGLScreenWidth,5.0*openGLScreenHeight/6.0);
				glTexCoord2f(1.0,0.0);
				glVertex2f(openGLScreenWidth,openGLScreenHeight/6.0);
			}
		}
		else if(type == HORIZONTAL)
		{
			glBegin(GL_QUADS);
			{
				glTexCoord2f(0.0,0.0);
				glVertex2f(openGLScreenWidth/6.0,openGLScreenHeight/2.0);
				glTexCoord2f(0.0,1.0);
				glVertex2f(openGLScreenWidth/6.0,openGLScreenHeight);
				glTexCoord2f(1.0,1.0);
				glVertex2f(5.0*openGLScreenWidth/6.0,openGLScreenHeight);
				glTexCoord2f(1.0,0.0);
				glVertex2f(5.0*openGLScreenWidth/6.0,openGLScreenHeight/2.0);
			}
		}
		
		glEnd(); 
		glDisable( GL_TEXTURE_2D );
		glDisable(GL_COLOR_MATERIAL);
	}
	
	//feature points
	if(currentFeatureInfo != NULL)
	{
		glPointSize(2.0);
		glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
		glBegin(GL_POINTS);
		for(int i=0; i<currentFeatureNo; i++)
		{
			Vec2f c(currentFeatureInfo[i].x, currentFeatureInfo[i].y);
			c /= scaleOfImg;
			c[0] = sxmin + ((sxmax-sxmin)*(c[0]-currentImageWidth*txmin)/(currentImageWidth*(txmax-txmin)));
			c[1] = symin + ((symax-symin)*(c[1]-currentImageHeight*tymin)/(currentImageHeight*(tymax-tymin)));
			if(c[0]>=sxmin && c[0]<=sxmax && c[1]>=symin && c[1]<=symax)
				glVertex2f(c[0],c[1]);
		}
		glEnd();
	}
	
	//projected points
	if(noProjectedPoints > 0)
	{
		glPointSize(2.0);
		glBegin(GL_POINTS);
		
		for(int i=0; i<noProjectedPoints; i++)
		{
			glColor4f(pointProjectionColors[i][0], pointProjectionColors[i][1], pointProjectionColors[i][2], 1.0f);
			
			Vec2f c(pointProjections[i*2], pointProjections[i*2+1]);
			c /= scaleOfImg;
			c[0] = sxmin + ((sxmax-sxmin)*(c[0]-currentImageWidth*txmin)/(currentImageWidth*(txmax-txmin)));
			c[1] = symin + ((symax-symin)*(c[1]-currentImageHeight*tymin)/(currentImageHeight*(tymax-tymin)));
			if(c[0]>=sxmin && c[0]<=sxmax && c[1]>=symin && c[1]<=symax)
				glVertex2f(c[0],c[1]);
			
		}
		glEnd();
	}
	
	//repetitions
	glEnable(GL_BLEND);
	glPointSize(5.0);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	for(int i=0; i<currentRepetitions.size(); i++)
	{
		//glColor4f(currentRepetitionColors[i][0], currentRepetitionColors[i][1], currentRepetitionColors[i][2], 1.0f);
		glColor4f(1.0, 0.0, 0.0, 0.5);
		glBegin(GL_POLYGON);
		for(int j=0; j<currentRepetitions[i].size(); j++)
		{
			Vec2f pos = currentRepetitions[i][j];
			pos /= scaleOfImg;
			pos[0] = sxmin + ((sxmax-sxmin)*(pos[0]-currentImageWidth*txmin)/(currentImageWidth*(txmax-txmin)));
			pos[1] = symin + ((symax-symin)*(pos[1]-currentImageHeight*tymin)/(currentImageHeight*(tymax-tymin)));
			if(pos[0]>=sxmin && pos[0]<=sxmax && pos[1]>=symin && pos[1]<=symax)
				glVertex2f(pos[0],pos[1]);
		}
		glEnd();
	}
	glDisable(GL_BLEND);
	
	glEnable(GL_BLEND);
	glPointSize(5.0);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	for(int i=0; i<currentRefinedRepetitions.size(); i++)
	{
		glColor4f(0.0, 0.0, 1.0, 0.5);
		glBegin(GL_POLYGON);
		for(int j=0; j<currentRefinedRepetitions[i].size(); j++)
		{
			Vec2f pos = currentRefinedRepetitions[i][j];
			pos /= scaleOfImg;
			pos[0] = sxmin + ((sxmax-sxmin)*(pos[0]-currentImageWidth*txmin)/(currentImageWidth*(txmax-txmin)));
			pos[1] = symin + ((symax-symin)*(pos[1]-currentImageHeight*tymin)/(currentImageHeight*(tymax-tymin)));
			if(pos[0]>=sxmin && pos[0]<=sxmax && pos[1]>=symin && pos[1]<=symax)
				glVertex2f(pos[0],pos[1]);
		}
		glEnd();
	}
	glDisable(GL_BLEND);
	
	//feature correspondences
	if(currentFeature != Vec2f(0.0, 0.0) && linkedFeature != Vec2f(0.0, 0.0))
	{
		Vec2f c,l;
		if(type == VERTICAL)
		{
			c[0] = currentFeature[0]*(openGLScreenWidth/2.0)/currentImageWidth;
			c[1] = currentFeature[1]*(4.0*openGLScreenHeight/6.0)/currentImageHeight + (openGLScreenHeight/6.0);
			l[0] = linkedFeature[0]*(openGLScreenWidth/2.0)/linkedImageWidth+openGLScreenWidth/2.0;
			l[1] = linkedFeature[1]*(4.0*openGLScreenHeight/6.0)/linkedImageHeight + (openGLScreenHeight/6.0);
		}
		else if(type == HORIZONTAL)
		{
			c[0] = currentFeature[0]*(4.0*openGLScreenWidth/6.0)/currentImageWidth + (openGLScreenWidth/6.0);
			c[1] = currentFeature[1]*(openGLScreenHeight/2.0)/currentImageHeight;
			l[0] = linkedFeature[0]*(4.0*openGLScreenWidth/6.0)/linkedImageWidth + (openGLScreenWidth/6.0);
			l[1] = linkedFeature[1]*(openGLScreenHeight/2.0)/linkedImageHeight + (openGLScreenHeight/2.0);
		}
	
		printf("c:%f %f, l:%f %f\n", c[0], c[1], l[0], l[1]);
	
		glLineWidth(1.0);
		glEnable(GL_LINE_SMOOTH);
		glBegin(GL_LINES);
	
		if(isValidCorr)
			glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		else
			glColor4f(1.0f, 0.0f, 0.0f, 0.0f);
		glVertex2f(c[0],c[1]);
		glVertex2i(l[0],l[1]);
	
		glEnd();
		glDisable(GL_LINE_SMOOTH);
	
		glPointSize(3.0);
		glBegin(GL_POINTS);
		glColor4f(1.0f, 0.0f, 0.0f, 0.0f);
		glVertex2f(c[0],c[1]);
		glVertex2i(l[0],l[1]);
		glEnd();
		
		if(epipolarLine.size() == 2)
		{
			Vec2f s,e;
			if(type == VERTICAL)
			{
				s[0] = epipolarLine[0][0]*(openGLScreenWidth/2.0)/linkedImageWidth+openGLScreenWidth/2.0;
				s[1] = epipolarLine[0][1]*(4.0*openGLScreenHeight/6.0)/linkedImageHeight + (openGLScreenHeight/6.0);
				e[0] = epipolarLine[1][0]*(openGLScreenWidth/2.0)/linkedImageWidth+openGLScreenWidth/2.0;
				e[1] = epipolarLine[1][1]*(4.0*openGLScreenHeight/6.0)/linkedImageHeight + (openGLScreenHeight/6.0);
			}
			else if(type == HORIZONTAL)
			{
				c[0] = epipolarLine[0][0]*(4.0*openGLScreenWidth/6.0)/linkedImageWidth + (openGLScreenWidth/6.0);
				c[1] = epipolarLine[0][1]*(openGLScreenHeight/2.0)/linkedImageHeight;
				l[0] = epipolarLine[1][0]*(4.0*openGLScreenWidth/6.0)/linkedImageWidth + (openGLScreenWidth/6.0);
				l[1] = epipolarLine[1][1]*(openGLScreenHeight/2.0)/linkedImageHeight + (openGLScreenHeight/2.0);
			}
			
			glLineWidth(1.0);
			glEnable(GL_LINE_SMOOTH);
			glBegin(GL_LINES);
			
			glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
			glVertex2f(s[0],s[1]);
			glVertex2i(e[0],e[1]);
			
			glEnd();
			glDisable(GL_LINE_SMOOTH);
		}
	
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
	
	}
	if(currentFeatures.size() > 0 && linkedFeatures.size() > 0)
	{
		Vec2f c,l;
		glLineWidth(1.0);
		glEnable(GL_LINE_SMOOTH);
		
		if(type == VERTICAL)
		{
			for(int i=0; i<currentFeatures.size(); i++)
			{
				c[0] = currentFeatures[i][0]*(openGLScreenWidth/2.0)/currentImageWidth;
				c[1] = currentFeatures[i][1]*(4.0*openGLScreenHeight/6.0)/currentImageHeight + (openGLScreenHeight/6.0);
				l[0] = linkedFeatures[i][0]*(openGLScreenWidth/2.0)/linkedImageWidth+openGLScreenWidth/2.0;
				l[1] = linkedFeatures[i][1]*(4.0*openGLScreenHeight/6.0)/linkedImageHeight + (openGLScreenHeight/6.0);
				
				if(validCorr[i])
					glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
				else
					glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
				
				glBegin(GL_LINES);
				
				glVertex2f(c[0],c[1]);
				glVertex2i(l[0],l[1]);
				
				glEnd();
			}
		}
		else if(type == HORIZONTAL)
		{
			for(int i=0; i<currentFeatures.size(); i++)
			{
				c[0] = currentFeatures[i][0]*(4.0*openGLScreenWidth/6.0)/currentImageWidth + (openGLScreenWidth/6.0);
				c[1] = currentFeatures[i][1]*(openGLScreenHeight/2.0)/currentImageHeight;
				l[0] = linkedFeatures[i][0]*(4.0*openGLScreenWidth/6.0)/linkedImageWidth + (openGLScreenWidth/6.0);
				l[1] = linkedFeatures[i][1]*(openGLScreenHeight/2.0)/linkedImageHeight + (openGLScreenHeight/2.0);
				
				if(validCorr[i])
					glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
				else
					glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
				
				glBegin(GL_LINES);
				
				glVertex2f(c[0],c[1]);
				glVertex2i(l[0],l[1]);
				
				glEnd();
			}
		}
			
		glDisable(GL_LINE_SMOOTH);
		
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	
	glFinish();
}

- (void)drawRect:(NSRect)dirtyRect {
    [self display];
}

-(void)makeCurrent
{
	NSOpenGLContext * glcontext;
	glcontext = [self openGLContext];
	[glcontext makeCurrentContext];
}

-(void)reshape
{
    glGetIntegerv (GL_VIEWPORT, viewportDimensions);
	
	openGLScreenWidth = viewportDimensions[2];
	openGLScreenHeight = viewportDimensions[3];
	
	glViewport(0, 0, openGLScreenWidth, openGLScreenHeight);
	
	[self makeCurrent];
	
	glGetIntegerv (GL_VIEWPORT, viewportDimensions);
}


-(void) attachBundlerParser:(BundlerOutputParser *)_bundleParser
{
	bundleParser = _bundleParser;
}

-(void) addNewImagePath:(std::string) filename
{
	listOfInputImagePaths.push_back(filename);
}

- (void) setScaleOfImg:(float)_scale
{
	scaleOfImg = _scale;
}

- (void) setFeaturePointInfo:(keypt_t*)featurePointInfo 
			   withFeatureNo:(int)noOfFeatures
{
	currentFeatureInfo = featurePointInfo;
	currentFeatureNo = noOfFeatures;
	currentFeature = Vec2f(0.0, 0.0);
	linkedFeature = Vec2f(0.0, 0.0);
	epipolarLine.clear();
}

-(void) setProjectedPointInfo:(float*)projectedPoints
					withColor:(Color*)projectedPointColors
				  withPointNo:(int)noOfPoints
{
	pointProjections = projectedPoints;
	pointProjectionColors = projectedPointColors;
	noProjectedPoints = noOfPoints;
}

- (void) setRepetitions:(vector<vector<Vec2f> > )_repetitions
			 withColors:(vector<Color>)_colors
{
	for(int i=0; i<currentRepetitions.size(); i++)
		currentRepetitions[i].clear();
	currentRepetitions.clear();
	currentRepetitionColors.clear();
	currentRepetitions = _repetitions;
	currentRepetitionColors = _colors;
}

-(void) setRefinedRepetitions:(vector<vector<Vec2f> >)_refinedRepetitions
{
	for(int i=0; i<currentRefinedRepetitions.size(); i++)
		currentRefinedRepetitions[i].clear();
	currentRefinedRepetitions.clear();
	currentRefinedRepetitions = _refinedRepetitions;
}

- (void) setGridStrokes: (vector<vector<Vec2f> >) _gridStrokes
{
	currentFeatures.clear();
	
	for(int i=0; i<currentGridStrokes.size(); i++)
		currentGridStrokes[i].clear();
	currentGridStrokes.clear();
	currentGridStrokes = _gridStrokes;
}

-(void) updateActiveImage: (int) _currentImageIndex
{
	showRectifiedImage = false;
	drawLine = false;
	lineScreenCoordinates.clear();
	noProjectedPoints = 0;
	for(int i=0; i<currentRepetitions.size(); i++)
		currentRepetitions[i].clear();
	currentRepetitions.clear();
	if(linkedImageBitmap!=NULL)
	{
		delete [] linkedImageBitmap;
		linkedImageBitmap = NULL;
	}
	
	imgZoomScale = 1.0;
	imgTranslate = Vec2f(0.0, 0.0);
	maxImgTranslate = (1.0 - imgZoomScale)/2.0;
	imgTranslateDelta = maxImgTranslate / 5.0;
	
	currentImageIndex = _currentImageIndex;
	printf("currentImage:%d\n", currentImageIndex);
	if(allImagesLoaded)
	{
		currentImageBitmap = allInputImages[currentImageIndex];
	}
	else
	{
		if(currentImageBitmap!=NULL)
		{
			delete [] currentImageBitmap;
			//free(currentImageBitmap);
			currentImageBitmap = NULL;
		}
		[self loadImage:currentImageIndex withWidth:currentImageWidth withHeight:currentImageHeight withData:&currentImageBitmap withScale:scaleOfImg];
	}
}

-(void) updateRectifiedImage: (Img*) _rectifiedImage
{
	showRectifiedImage = true;
	drawLine = false;
	lineScreenCoordinates.clear();
	noProjectedPoints = 0;
	for(int i=0; i<currentRepetitions.size(); i++)
		currentRepetitions[i].clear();
	currentRepetitions.clear();
	
	imgZoomScale = 1.0;
	imgTranslate = Vec2f(0.0, 0.0);
	maxImgTranslate = (1.0 - imgZoomScale)/2.0;
	imgTranslateDelta = maxImgTranslate / 5.0;
	
	rectImageWidth = _rectifiedImage->width();
	rectImageHeight = _rectifiedImage->height();
	if(rectifiedImage != NULL)
		delete [] rectifiedImage;
	
	rectifiedImage = new float[rectImageWidth*rectImageHeight*4];
	
	for(int y=0; y<rectImageHeight; y++)
	{
		for(int x=0; x<rectImageWidth; x++)
		{
			*(rectifiedImage + 4*((rectImageHeight-y-1) * rectImageWidth + x)) = (*_rectifiedImage)(x,rectImageHeight-y-1)[0];		// red
			*(rectifiedImage + 4*((rectImageHeight-y-1) * rectImageWidth + x)+1) = (*_rectifiedImage)(x,rectImageHeight-y-1)[1];	// green
			*(rectifiedImage + 4*((rectImageHeight-y-1) * rectImageWidth + x)+2) = (*_rectifiedImage)(x,rectImageHeight-y-1)[2];	// blue
			*(rectifiedImage + 4*((rectImageHeight-y-1) * rectImageWidth + x)+3) = 255;	// alpha
		}
	}
	
	rectImageScale = 1.0;
	while(rectImageWidth/rectImageScale > openGLScreenWidth || rectImageHeight/rectImageScale > openGLScreenHeight)
	{
		rectImageScale *= 2.0;
	}
}

-(void) updateActiveImages: (int) _currentImageIndex
				  linkedTo: (int) _linkedImageIndex
{
	showRectifiedImage = false;
	drawLine = false;
	lineScreenCoordinates.clear();
	noProjectedPoints = 0;
	for(int i=0; i<currentRepetitions.size(); i++)
		currentRepetitions[i].clear();
	currentRepetitions.clear();
	
	imgZoomScale = 1.0;
	imgTranslate = Vec2f(0.0, 0.0);
	maxImgTranslate = (1.0 - imgZoomScale)/2.0;
	imgTranslateDelta = maxImgTranslate / 5.0;
	
	currentImageIndex = _currentImageIndex;
	linkedImageIndex = _linkedImageIndex;
	
	printf("currentImage:%d, linkedImage:%d\n", currentImageIndex, linkedImageIndex);
	
	if(allImagesLoaded)
	{
		currentImageBitmap = allInputImages[currentImageIndex];
		linkedImageBitmap = allInputImages[linkedImageIndex];
	}
	
	else
	{
		if(currentImageBitmap!=NULL)
		{
			delete [] currentImageBitmap;
			//free(currentImageBitmap);
			currentImageBitmap = NULL;
		}
		[self loadImage:currentImageIndex withWidth:currentImageWidth withHeight:currentImageHeight withData:&currentImageBitmap withScale:scaleOfImg];
	
		if(linkedImageBitmap!=NULL)
		{
			delete [] linkedImageBitmap;
			linkedImageBitmap = NULL;
		}
		[self loadImage:linkedImageIndex withWidth:linkedImageWidth withHeight:linkedImageHeight withData:&linkedImageBitmap withScale:scaleOfImg];
	}
}

-(void) updateActiveRectImages: (Img *) _currentImage 
					  linkedTo: (Img *) _linkedImage
{
	showRectifiedImage = false;
	drawLine = false;
	lineScreenCoordinates.clear();
	noProjectedPoints = 0;
	for(int i=0; i<currentRepetitions.size(); i++)
		currentRepetitions[i].clear();
	currentRepetitions.clear();
	
	imgZoomScale = 1.0;
	imgTranslate = Vec2f(0.0, 0.0);
	maxImgTranslate = (1.0 - imgZoomScale)/2.0;
	imgTranslateDelta = maxImgTranslate / 5.0;
	
	if(currentImageBitmap!=NULL)
	{
		delete [] currentImageBitmap;
		currentImageBitmap = NULL;
	}
	
	int w = _currentImage->width();
	int h = _currentImage->height();
	
	currentImageBitmap = new float[w*h*4];
	
	for (int y = 0; y < h; ++y) 
	{      
		for (int x = 0; x < w; ++x) 
		{
			*(currentImageBitmap + 4*((h-y-1) * w + x)) = (*_currentImage)(x, h-y-1)[0];
			*(currentImageBitmap + 4*((h-y-1) * w + x)+1) = (*_currentImage)(x, h-y-1)[1];
			*(currentImageBitmap + 4*((h-y-1) * w + x)+2) = (*_currentImage)(x, h-y-1)[2];
			*(currentImageBitmap + 4*((h-y-1) * w + x)+3) = 1.0;
		}
	}
	
	if(linkedImageBitmap!=NULL)
	{
		delete [] linkedImageBitmap;
		linkedImageBitmap = NULL;
	}
	
	w = _linkedImage->width();
	h = _linkedImage->height();
	
	linkedImageBitmap = new float[w*h*4];
	
	for (int y = 0; y < h; ++y) 
	{      
		for (int x = 0; x < w; ++x) 
		{
			*(linkedImageBitmap + 4*((h-y-1) * w + x)) = (*_linkedImage)(x, h-y-1)[0];
			*(linkedImageBitmap + 4*((h-y-1) * w + x)+1) = (*_linkedImage)(x, h-y-1)[1];
			*(linkedImageBitmap + 4*((h-y-1) * w + x)+2) = (*_linkedImage)(x, h-y-1)[2];
			*(linkedImageBitmap + 4*((h-y-1) * w + x)+3) = 1.0;
		}
	}
	
}

-(void)updateActiveCorrespondences:(Vec2f) _currentFeature
						  linkedTo:(Vec2f) _linkedFeature
				 withEpipolarStart:(Vec2f) _start
				   withEpipolarEnd:(Vec2f) _end
{
	noProjectedPoints = 0;
	currentFeature = _currentFeature;
	linkedFeature = _linkedFeature;
	epipolarLine.clear();
	epipolarLine.push_back(_start);
	epipolarLine.push_back(_end);
	currentFeatureInfo = NULL;
	currentFeatureNo = 0;
	printf("currentFeature:%f %f, linkedFeature:%f %f\n", currentFeature[0], currentFeature[1], linkedFeature[0], linkedFeature[1]);
}

-(void)updateCorrespondences:(vector<Vec2f>) _currentFeatures
					linkedTo:(vector<Vec2f>) _linkedFeatures
					validity:(vector<bool>) _validCorr
{
	noProjectedPoints = 0;
	currentFeatures.clear();
	linkedFeatures.clear();
	validCorr.clear();
	
	currentFeatures = _currentFeatures;
	linkedFeatures = _linkedFeatures;
	validCorr = _validCorr;
	currentFeatureInfo = NULL;
	currentFeatureNo = 0;
}

-(void)updateActiveCorrespondences:(Vec2f) _currentFeature
						  linkedTo:(Vec2f) _linkedFeature
						   isValid:(bool) _valid
{
	noProjectedPoints = 0;
	currentFeature = _currentFeature;
	linkedFeature = _linkedFeature;
	isValidCorr = _valid;
	epipolarLine.clear();
	currentFeatureInfo = NULL;
	currentFeatureNo = 0;
	printf("currentFeature:%f %f, linkedFeature:%f %f\n", currentFeature[0], currentFeature[1], linkedFeature[0], linkedFeature[1]);
}

-(void) loadImage: (int) imageIndex
		withWidth:(int &) w
	   withHeight:(int &) h
		 withData:(float**) data
		withScale:(float) scale
{
	Img currentImg;
	currentImg.read(listOfInputImagePaths[imageIndex].c_str());
	w = currentImg.width();
	h = currentImg.height();
	
	//apply gaussian filter
	std::vector<float> mask;
	mask.resize(16, 0.0);
	mask[0] = 1.0/64.0; mask[1] = 3.0/64.0; mask[2] = 3.0/64.0; mask[3] = 1.0/64.0;
	mask[4] = 3.0/64.0; mask[5] = 9.0/64.0; mask[6] = 9.0/64.0; mask[7] = 3.0/64.0;
	mask[8] = 3.0/64.0; mask[9] = 9.0/64.0; mask[10] = 9.0/64.0; mask[11] = 3.0/64.0;
	mask[12] = 1.0/64.0; mask[13] = 3.0/64.0; mask[14] = 9.0/64.0; mask[15] = 1.0/64.0;
	
	int wNew = w / scale;
	int hNew = h / scale;
	(*data) = new float[wNew*hNew*4];
	
	for (int y = 0; y < hNew; ++y) 
	{      
		for (int x = 0; x < wNew; ++x) 
		{
			if(scale != 1.0)
			{
				float rTotal = 0.0;
				float gTotal = 0.0;
				float bTotal = 0.0;
				float denom = 0.0;
				for (int j = -1; j < 3; ++j) 
				{
					int ytmp = (int)scale * y + j;
					if (ytmp < 0 || h - 1 < ytmp)
						continue;
					
					for (int i = -1; i < 3; ++i) 
					{
						int xtmp = (int)scale * x + i;
						if (xtmp < 0 || w - 1 < xtmp)
							continue;
						
						float r = currentImg(xtmp, h-ytmp-1)[0];
						float g = currentImg(xtmp, h-ytmp-1)[1];
						float b = currentImg(xtmp, h-ytmp-1)[2];
						
						r = mask[i+1 + (j+1)*4] * r;
						g = mask[i+1 + (j+1)*4] * g;
						b = mask[i+1 + (j+1)*4] * b;
						
						rTotal += r;
						gTotal += g;
						bTotal += b;
						
						denom += mask[i+1 + (j+1)*4];
					}
					
				}
				rTotal /= denom; gTotal /= denom; bTotal /= denom;
				*((*data) + 4*((hNew-y-1) * wNew + x)) = rTotal;		// red
				*((*data) + 4*((hNew-y-1) * wNew + x)+1) = gTotal;	// green
				*((*data) + 4*((hNew-y-1) * wNew + x)+2) = bTotal;	// blue
				*((*data) + 4*((hNew-y-1) * wNew + x)+3) = 255;	// alpha
			}
			else
			{
				*((*data) + 4*((hNew-y-1) * wNew + x)) = currentImg(x, hNew-y-1)[0];
				*((*data) + 4*((hNew-y-1) * wNew + x)+1) = currentImg(x, hNew-y-1)[1];
				*((*data) + 4*((hNew-y-1) * wNew + x)+2) = currentImg(x, hNew-y-1)[2];
				*((*data) + 4*((hNew-y-1) * wNew + x)+3) = 1.0;
			}
		}
	}

	w = wNew; h = hNew;
	
}

-(void) loadAllImages:(float) scale
{
	for(int i=0; i<allInputImages.size(); i++)
	{
		if(allInputImages[i] != NULL)
		{
			delete [] allInputImages[i];
			allInputImages[i] = NULL;
		}
	}
	allInputImages.clear();
	
	int w, h;
	for(int i=0; i<listOfInputImagePaths.size(); i++)
	{
		allInputImages.push_back(new float());
		[self loadImage: i withWidth:w withHeight:h withData:&(allInputImages[i]) withScale:scale];
	}
	currentImageWidth = w; currentImageHeight = h;
	linkedImageWidth = w; linkedImageHeight = h;
	
	allImagesLoaded = true;
}

- (void) takeScreenshot
{
	std::string newFileName("imageOfFrame");
	
	char counterBuffer[16];
	snprintf(counterBuffer,16,"%06d",numberOfScreenshotsTaken);
	newFileName.append(counterBuffer);
	newFileName.append(std::string(".png"));
	
	NSString * newNSFileName = [NSString stringWithUTF8String: newFileName.c_str()];
	[self takeScreenshot:newNSFileName];
	numberOfScreenshotsTaken++;
}

-(void)takeScreenshot:(NSString *) _fileName
{
	NSAutoreleasePool * pool = [[NSAutoreleasePool alloc] init];
	[self makeCurrent];
	
	int imageWidth	= viewportDimensions[2];
	int imageHeight = viewportDimensions[3];
	
	NSLog(@"width: %d, height: %d",imageWidth,imageHeight);
	
	float * data = new float[imageWidth * imageHeight*4]; 	
	
	[self renderOffScreenWithWidth:imageWidth andHeight: imageHeight andData: data];
	
	NSImage * screenshotImage;
	NSBitmapImageRep * imageRepresentation;
	NSData * imageData;
	
	// Create image
	screenshotImage=[[NSImage alloc] initWithSize:NSMakeSize(imageWidth, imageHeight)];
	
	imageRepresentation = [[NSBitmapImageRep alloc]
						   initWithBitmapDataPlanes:NULL
						   pixelsWide:imageWidth 
						   pixelsHigh:imageHeight
						   bitsPerSample:8 
						   samplesPerPixel:4
						   hasAlpha:YES
						   isPlanar:NO
						   colorSpaceName:NSCalibratedRGBColorSpace
						   bytesPerRow:imageWidth*4 // assuming no empty space at the end of a row 
						   bitsPerPixel:NULL // without meaningless bits
						   ]; 
	
	unsigned char * imageBitmapData=[imageRepresentation bitmapData];
	
	
	// copy data to imageBitmapData with vertical inversion
	
	unsigned int i=0;
	for(unsigned int y=0;y<imageHeight;y++)
	{
		for(unsigned int x=0;x<imageWidth;x++,i+=4)
		{
			*(imageBitmapData + 4* ((imageHeight-y-1) * imageWidth + x))   = data[i]  *255.0;	// red
			*(imageBitmapData + 4* ((imageHeight-y-1) * imageWidth + x)+1) = data[i+1]*255.0;	// green
			*(imageBitmapData + 4* ((imageHeight-y-1) * imageWidth + x)+2) = data[i+2]*255.0;	// blue
			*(imageBitmapData + 4* ((imageHeight-y-1) * imageWidth + x)+3) = data[i+3]*255.0;	// alpha
		}
		
	}
	
	[screenshotImage addRepresentation:imageRepresentation];	
	
	// use file manager to check if Screenshots directory exists
	
	NSFileManager * newFileManager = [NSFileManager defaultManager];
	
	if ([newFileManager fileExistsAtPath:[[NSHomeDirectory() stringByAppendingPathComponent:@"Desktop"] stringByAppendingString:@"/Screenshots/"]]== NO)
	{
		[newFileManager createDirectoryAtPath:[[NSHomeDirectory() stringByAppendingPathComponent:@"Desktop"] stringByAppendingString:@"/Screenshots/"] attributes:nil];
	}
	
	imageData = [imageRepresentation representationUsingType: NSPNGFileType properties: nil];
	
	[imageData writeToFile: [[[NSHomeDirectory() stringByAppendingPathComponent:@"Desktop"] stringByAppendingString:@"/Screenshots/"] stringByAppendingString:_fileName] atomically: NO];
	
	//delete imageBitmapData;
	
	[imageRepresentation release];
	[screenshotImage release];	
	delete [] data;
	
	[pool release];
}

-(void) renderOffScreenWithWidth:(unsigned int) width andHeight: (unsigned int) height andData: (float *) data
{
	if(data!=NULL)
	{
		GLuint renderbuffer;
		GLenum status;
		
		// Set up a FBO with one renderbuffer attachment
		glGenRenderbuffersEXT(1, &renderbuffer);
		
		glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, renderbuffer);
		
		glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA, width, height);
		
		glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
									 GL_RENDERBUFFER_EXT, renderbuffer);
		
		status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
		
		if (status != GL_FRAMEBUFFER_COMPLETE_EXT)
		{
			// Handle errors
			printf("Framebuffer is incomplete\n");
		}
		
		//Your code to draw content to the renderbuffer
		// ...
		
		float alphaBackup = backgroundColor[3];
		backgroundColor[3] = 0.0;
		[self display];
		backgroundColor[3] = alphaBackup;
		
		// Make the window the target
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		
		//Your code to use the contents
		
		glReadBuffer(GL_COLOR_ATTACHMENT0_EXT); 
		glReadPixels(0, 0, width, height, GL_RGBA, GL_FLOAT, data); 
		
		// Delete the renderbuffer attachment
		glDeleteRenderbuffersEXT(1, &renderbuffer);
	}
}

-(void) setGrids3DComputed:(bool)_grids3DComputed
{
	grids3DComputed = _grids3DComputed;
}

-(void) setDrawLine:(bool)_drawLine
{
	if(showRectifiedImage)
		drawLine = _drawLine;
	else if(grids3DComputed)
		drawLine = _drawLine;
	
	lineScreenCoordinates.clear();
}

-(void) translateImageLeft
{
	imgTranslate[0] -= imgTranslateDelta;
	if(imgTranslate[0] < -maxImgTranslate)
		imgTranslate[0] = -maxImgTranslate;
}

-(void) translateImageRight
{
	imgTranslate[0] += imgTranslateDelta;
	if(imgTranslate[0] > maxImgTranslate)
		imgTranslate[0] = maxImgTranslate;
}

-(void) translateImageUp
{
	imgTranslate[1] -= imgTranslateDelta;
	if(imgTranslate[1] < -maxImgTranslate)
		imgTranslate[1] = -maxImgTranslate;
}

-(void) translateImageDown
{
	imgTranslate[1] += imgTranslateDelta;
	if(imgTranslate[1] > maxImgTranslate)
		imgTranslate[1] = maxImgTranslate;
}

-(void) scaleImage
{
	imgZoomScale -= 0.2;
	if(imgZoomScale <= 0.0)
		imgZoomScale = 1.0;
	imgTranslate[0] = 0.0;
	imgTranslate[1] = 0.0;
	maxImgTranslate = (1.0-imgZoomScale)/2.0;
	imgTranslateDelta = maxImgTranslate/5.0;
}

-(void) mouseDown:(NSEvent *)theEvent
{
	NSPoint pt = [self convertPoint:[theEvent locationInWindow] fromView:nil];
	printf("pt:%f %f\n", pt.x, pt.y);
	
	if(drawLine)
	{
		startLocation = pt;
		startLocation.y = openGLScreenHeight-startLocation.y;
		lineScreenCoordinates.push_back(Vec2f(startLocation.x, startLocation.y));
		lineScreenCoordinates.push_back(Vec2f(startLocation.x, startLocation.y));
	}
	else
	{
		float x = pt.x;
		float y = openGLScreenHeight - pt.y;
		if(x>=sxmin && x<=sxmax && y>=symin && y<=symax)
		{
			float ratioX = (x-sxmin)/(sxmax-sxmin);
			float ratioY = (y-symin)/(symax-symin);
			float tX = txmin + ratioX*(txmax-txmin);
			float tY = tymin + ratioY*(tymax-tymin);
			currentCoord[0] = currentImageWidth*tX;
			currentCoord[1] = currentImageHeight*tY;
			printf("coord:%f %f\n", currentCoord[0], currentCoord[1]);
		}
	}
}

- (void)mouseDragged:(NSEvent *)theEvent
{
	if(drawLine)
	{
		NSPoint mouseLocation = [self convertPoint:[theEvent locationInWindow] fromView:nil];
		mouseLocation.y = openGLScreenHeight-mouseLocation.y;
		lineScreenCoordinates[lineScreenCoordinates.size()-1][0] = mouseLocation.x;
		lineScreenCoordinates[lineScreenCoordinates.size()-1][1] = mouseLocation.y;
		
		[self setNeedsDisplay:YES];
	}
}

- (void)mouseUp:(NSEvent *)theEvent
{
	if(drawLine)
	{
		NSPoint mouseLocation = [self convertPoint:[theEvent locationInWindow] fromView:nil];
		mouseLocation.y = openGLScreenHeight-mouseLocation.y;
		lineScreenCoordinates[lineScreenCoordinates.size()-1][0] = mouseLocation.x;
		lineScreenCoordinates[lineScreenCoordinates.size()-1][1] = mouseLocation.y;
		[self setNeedsDisplay:YES];
	}
}

-(void) getTemplateCoord:(Vec2f&) topLeft
				   width:(int&)w
				  height:(int&)h
{
	float xmin, xmax, ymin, ymax;
	topLeft[0] = 0.0; topLeft[1] = 0.0;
	w = 0;
	h = 0;
	int noOfPaths = lineScreenCoordinates.size();
	
	for(int i=0; i<noOfPaths; i++)
	{
		if(i==0)
		{
			xmin = lineScreenCoordinates[0][0];
			xmax = lineScreenCoordinates[0][0];
			ymin = lineScreenCoordinates[0][1];
			ymax = lineScreenCoordinates[0][1];
		}
		if(lineScreenCoordinates[i][0] < xmin)
			xmin = lineScreenCoordinates[i][0];
		else if(lineScreenCoordinates[i][0] > xmax)
			xmax = lineScreenCoordinates[i][0];
		
		if(lineScreenCoordinates[i][1] < ymin)
			ymin = lineScreenCoordinates[i][1];
		else if(lineScreenCoordinates[i][1] > ymax)
			ymax = lineScreenCoordinates[i][1];
	}
	
	Vec2f start((openGLScreenWidth - rectImageWidth/rectImageScale)/2.0, (openGLScreenHeight-rectImageHeight/rectImageScale)/2.0);
	Vec2f end(start[0]+rectImageWidth/rectImageScale, start[1]+rectImageHeight/rectImageScale);
	
	float ratioX = (xmin-start[0])/(end[0]-start[0]);
	float ratioY = (ymin-start[1])/(end[1]-start[1]);
	float tX = txmin + ratioX*(txmax-txmin);
	float tY = tymin + ratioY*(tymax-tymin);
	topLeft[0] = rectImageWidth*tX;
	topLeft[1] = rectImageHeight*tY;
	
	ratioX = (xmax-start[0])/(end[0]-start[0]);
	ratioY = (ymax-start[1])/(end[1]-start[1]);
	tX = txmin + ratioX*(txmax-txmin);
	tY = tymin + ratioY*(tymax-tymin);
	
	w = floor(rectImageWidth*tX - topLeft[0] + 0.5);
	h = floor(rectImageHeight*tY - topLeft[1] + 0.5);
	
}

-(void) getStrokes:(vector<Vec2f>&) endpoints
{
	/*lineScreenCoordinates.push_back(Vec2f(541.000000, openGLScreenHeight-537.000000));
	lineScreenCoordinates.push_back(Vec2f(542.000000, openGLScreenHeight-505.000000));
	lineScreenCoordinates.push_back(Vec2f(542.000000, openGLScreenHeight-505.000000));
	lineScreenCoordinates.push_back(Vec2f(598.000000, openGLScreenHeight-514.000000));
	lineScreenCoordinates.push_back(Vec2f(598.000000, openGLScreenHeight-514.000000));
	lineScreenCoordinates.push_back(Vec2f(595.000000, openGLScreenHeight-547.000000));
	lineScreenCoordinates.push_back(Vec2f(595.000000, openGLScreenHeight-547.000000));
	lineScreenCoordinates.push_back(Vec2f(541.000000, openGLScreenHeight-537.000000));*/
	
	int noOfPaths = lineScreenCoordinates.size();
	
	float ratioX = currentImageWidth * (txmax-txmin) / (sxmax-sxmin);
	float ratioY = currentImageHeight * (tymax-tymin) / (symax-symin);
	
	for(int i=0; i<noOfPaths; i++)
	{
		if(i%2 == 0)
			continue;
		
		float x = lineScreenCoordinates[i][0];
		float y = lineScreenCoordinates[i][1];
		x = (x-sxmin) *ratioX + currentImageWidth*txmin;
		x *= scaleOfImg;
		y = (y-symin) * ratioY + currentImageHeight*tymin;
		y *= scaleOfImg;
		endpoints.push_back(Vec2f(x,y));
	}
}

-(void) getCurrentCoord:(Vec2f &)c
{
	c = currentCoord;
	currentCoord[0] = 0.0;
	currentCoord[1] = 0.0;
}

@end
