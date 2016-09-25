/*
 *  SIFT.h
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 12/5/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SIFT_H
#define SIFT_H

#include<valarray>
#include<vector>
#include<ostream>
#include<cmath>
#include<limits>
#include "Common.h"
#include "Image.h"

typedef float pixel_t ;

struct PgmBuffer
{
	int width ;     ///< Image width
	int height ;    ///< Image hegith
	pixel_t* data ; ///< Image data
} ;

class SIFT
{
public:
	typedef std::vector<keypt_t>	Keypoints ;          ///< Keypoint list datatype
	typedef Keypoints::iterator       KeypointsIter ;      ///< Keypoint list iter datatype
	typedef Keypoints::const_iterator KeypointsConstIter ; ///< Keypoint list const iter datatype
		
	SIFT();
	SIFT(float _sigman, float _sigma0, int _O, int _S, int _omin, int _smin, int _smax) ;
	~SIFT() ;
	
	void prepareImage(Img *im);
	void cleanImage();
	
	void processImage(Img *im, std::vector<keypt_t> &keyInfo, std::vector<short> &keyDescr);
	void processUprightImage(Img *im, std::vector<keypt_t> &keyInfo, std::vector<short> &keyDescr);
	
	void extractPgm(Img *im, PgmBuffer& buffer);
	
	void process(const pixel_t* _im_pt, int _width, int _height) ;
		
	//octave querying
	float* getOctave(int o) ;
	float* getLevel(int o, int s) ;
	int		getWidth() const {return width; }
	int		getHeight() const {return height; }
	int		getOctaveWidth(int o) const ;
	int		getOctaveHeight(int o) const ;
	float	getOctaveSamplingPeriod(int o) const ;
	float	getScaleFromIndex(float o, float s) const ;
	keypt_t getKeypoint(float x,  float y, float s);
	
	//descriptor
	bool getNormalizeDescriptor() const {return normalizeDescriptor ;}
	void setNormalizeDescriptor(bool) ;
	void setMagnification(float) ;
	float getMagnification() const {return magnif ;}  
	
	//detector and descriptor
	void detectKeypoints(float threshold, float edgeThreshold) ;
	int computeKeypointOrientations(float angles [4], keypt_t keypoint) ; 
	void computeKeypointDescriptor(float* descrPt, keypt_t keypoint, float angle) ;
	KeypointsIter keypointsBegin() ;
	KeypointsIter keypointsEnd() ;
		
private:
	float   fast_resqrt(float x) ;
	double  fast_resqrt(double x) ;
	float fast_expn(float_t x) ;
	float fast_abs(float_t x) ;
	float fast_mod_2pi(float_t x) ;
	float fast_atan2(float_t y, float_t x) ;
	float fast_sqrt(float_t x) ;
	int	 fast_floor(float_t x) ;
	
	void copy(pixel_t* dst, pixel_t const* src, int width, int height);
	void copyAndUpsampleRows(pixel_t* dst, pixel_t const* src, int width, int height);
	void copyAndDownsample(pixel_t* dst, pixel_t const* src, int width, int height, int d);
	
	void normalize(pixel_t* filter, int W);
	void econvolve(pixel_t*       dst_pt, 
				   const pixel_t* src_pt,    int M, int N,
				   const pixel_t* filter_pt, int W);
	
	void prepareBuffers() ;
	void freeBuffers() ;
	void smooth(pixel_t       * dst, 
				pixel_t       * temp, 
				pixel_t const * src, 
				int width, int height, float_t s) ;
		
	void prepareGrad(int o) ;
	
	void normalizeHistogram(float* L_begin, float* L_end);
		
	// scale space parameters
	float sigman ;
	float sigma0 ;
	float sigmak ;
		
	int O ;
	int S ; 
	int omin ;
	int smin ; 
	int smax ;
		
	int width ;
	int height ;
		
	// descriptor parameters
	float magnif ;
	bool  normalizeDescriptor ;
		
	// buffers
	PgmBuffer	pgmBuffer;
	pixel_t*	temp ;
	int			tempReserved ;
	bool		tempIsGrad  ;
	int			tempOctave ;
	pixel_t**	octaves ;
	
	pixel_t*  filter ;
	int           filterReserved ;
		
	Keypoints keypoints ;  
} ;

#endif
