/*
 *  SIFT.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 12/5/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "SIFT.h"
#include<algorithm>
#include<iostream>
#include<sstream>
#include<cassert>

#define CHECK_NEIGHBORS(CMP,SGN)                    \
( v CMP ## = SGN 0.8 * threshold &&     \
v CMP *(pt + xo) &&                   \
v CMP *(pt - xo) &&                   \
v CMP *(pt + so) &&                   \
v CMP *(pt - so) &&                   \
v CMP *(pt + yo) &&                   \
v CMP *(pt - yo) &&                   \
\
v CMP *(pt + yo + xo) &&              \
v CMP *(pt + yo - xo) &&              \
v CMP *(pt - yo + xo) &&              \
v CMP *(pt - yo - xo) &&              \
\
v CMP *(pt + xo      + so) &&         \
v CMP *(pt - xo      + so) &&         \
v CMP *(pt + yo      + so) &&         \
v CMP *(pt - yo      + so) &&         \
v CMP *(pt + yo + xo + so) &&         \
v CMP *(pt + yo - xo + so) &&         \
v CMP *(pt - yo + xo + so) &&         \
v CMP *(pt - yo - xo + so) &&         \
\
v CMP *(pt + xo      - so) &&         \
v CMP *(pt - xo      - so) &&         \
v CMP *(pt + yo      - so) &&         \
v CMP *(pt - yo      - so) &&         \
v CMP *(pt + yo + xo - so) &&         \
v CMP *(pt + yo - xo - so) &&         \
v CMP *(pt - yo + xo - so) &&         \
v CMP *(pt - yo - xo - so) )

#define Aat(i,j)     (A[(i)+(j)*3]) 

#define atd(dbinx,dbiny,dbint) *(dpt + (dbint)*binto + (dbiny)*binyo + (dbinx)*binxo)

int SIFT::getOctaveWidth(int o) const
{
	if(omin > o || o >= omin+O)
	{
		printf("[SIFT] getOctaveWidth:invalid argument\n");
		return -1;
	}
	return (o >= 0) ? (width >> o) : (width << -o) ;
}

int SIFT::getOctaveHeight(int o) const
{
	if(omin > o || o >= omin+O)
	{
		printf("[SIFT] getOctaveHeight:invalid argument\n");
		return -1;
	}
	return (o >= 0) ? (height >> o) : (height << -o) ;
}

pixel_t * SIFT::getOctave(int o) 
{
	if(omin > o || o >= omin+O)
	{
		printf("[SIFT] getOctave:invalid argument\n");
		return NULL;
	}
	return octaves[o-omin] ;
}

pixel_t * SIFT::getLevel(int o, int s) 
{
	if(omin > o || o > omin+O)
	{
		printf("[SIFT] getLevel:invalid argument\n");
		return NULL;
	}
	if(smin > s || s > smax)
	{
		printf("[SIFT] getLevel:invalid argument\n");
		return NULL;
	}
	return octaves[o - omin] + getOctaveWidth(o)*getOctaveHeight(o) * (s-smin) ;
}

float SIFT::getOctaveSamplingPeriod(int o) const
{
	return (o >= 0) ? (1 << o) : 1.0f / (1 << -o) ;
}

float_t SIFT::getScaleFromIndex(float o, float s) const
{
	return sigma0 * powf( 2.0f, o + s / S ) ;
}

SIFT::KeypointsIter SIFT::keypointsBegin()
{
	return keypoints.begin() ;
}

SIFT::KeypointsIter SIFT::keypointsEnd()
{
	return keypoints.end() ;
}

void SIFT::setNormalizeDescriptor(bool flag)
{
	normalizeDescriptor = flag ;
}

void SIFT::setMagnification(float _magnif)
{
	magnif = _magnif ;
}

float SIFT::fast_expn(float x)
{
	return exp(-x) ;
}

float SIFT::fast_mod_2pi(float x)
{
	return (x>=0) ? std::fmod(x, float(2*M_PI)) : 2*M_PI + std::fmod(x, float(2*M_PI)) ;
}

int SIFT::fast_floor(float x)
{
	return int( std::floor(x) ) ;
}

float SIFT::fast_abs(float x)
{
	return std::fabs(x) ; 
}

float SIFT::fast_atan2(float y, float x)
{
	return std::atan2(y,x) ;
}

float SIFT::fast_resqrt(float x)
{
	return float(1.0) / std::sqrt(x) ;
}

double SIFT::fast_resqrt(double x)
{
	return double(1.0) / std::sqrt(x) ;
}

float SIFT::fast_sqrt(float x)
{
	return std::sqrt(x) ;
}

void SIFT::copy(pixel_t* dst, pixel_t const* src, int width, int height)
{
	memcpy(dst, src, sizeof(pixel_t)*width*height)  ;
}
		
void SIFT::copyAndUpsampleRows(pixel_t* dst, pixel_t const* src, int width, int height)
{
	for(int y = 0 ; y < height ; ++y) 
	{
		pixel_t b, a ;
		b = a = *src++ ;
		for(int x = 0 ; x < width-1 ; ++x) 
		{
			b = *src++ ;
			*dst = a ;         dst += height ;
			*dst = 0.5*(a+b) ; dst += height ;
			a = b ;
		}
		*dst = b ; dst += height ;
		*dst = b ; dst += height ;
		dst += 1 - width * 2 * height ;
	}  
}		

void SIFT::copyAndDownsample(pixel_t* dst, pixel_t const* src, int width, int height, int d)
{
	for(int y = 0 ; y < height ; y+=d)
	{
		pixel_t const * srcrowp = src + y * width ;    
		for(int x = 0 ; x < width - (d-1) ; x+=d) 
		{     
			*dst++ = *srcrowp ;
			srcrowp += d ;
		}
	}
}

void SIFT::normalize(pixel_t* filter, int W)
{
	pixel_t  acc  = 0 ; 
	pixel_t* iter = filter ;
	pixel_t* end  = filter + 2*W+1 ;
	while(iter != end) 
		acc += *iter++ ;
	
	iter = filter ;
	while(iter != end) 
		*iter++ /= acc ;
}

void SIFT::econvolve(pixel_t*       dst_pt, 
					 const pixel_t* src_pt,    int M, int N,
					 const pixel_t* filter_pt, int W)
{
	// convolve along columns, save transpose
	// image is M by N 
	// buffer is N by M 
	// filter is (2*W+1) by 1
	for(int j = 0 ; j < N ; ++j) 
	{
		for(int i = 0 ; i < M ; ++i) 
		{
			pixel_t   acc = 0.0 ;
			pixel_t const * g = filter_pt ;
			pixel_t const * start = src_pt + (i-W) ;
			pixel_t const * stop  ;
			pixel_t   x ;
			
			// beginning
			stop = src_pt + std::max(0, i-W) ;
			x    = *stop ;
			while( start <= stop ) 
			{ 
				acc += (*g++) * x ; start++ ; 
			}
			
			// middle
			stop =  src_pt + std::min(M-1, i+W) ;
			while( start <  stop ) 
				acc += (*g++) * (*start++) ;
			
			// end
			x  = *start ;
			stop = src_pt + (i+W) ;
			while( start <= stop ) 
			{
				acc += (*g++) * x ; start++ ; 
			} 
			
			// save 
			*dst_pt = acc ; 
			dst_pt += N ;
			
		}
		// next column
		src_pt += M ;
		dst_pt -= M*N - 1 ;
	}
}
		
void SIFT::smooth(pixel_t* dst, pixel_t* temp, pixel_t const* src, int width, int height, float s)
{
	// make sure a buffer large enough has been allocated to hold the filter
	int W = int( ceil( float(4.0) * s ) ) ;
	if( ! filter ) 
	{
		filterReserved = 0 ;
	}
		
	if( filterReserved < W ) 
	{
		filterReserved = W ;
		if( filter ) delete [] filter ;
		filter = new pixel_t [ 2* filterReserved + 1 ] ;
	}
		
	// pre-compute filter
	for(int j = 0 ; j < 2*W+1 ; ++j) 
		filter[j] = pixel_t(std::exp(float(-0.5 * (j-W) * (j-W) / (s*s) ))) ;
		
	// normalize to one
	normalize(filter, W) ;
		
	// convolve
	econvolve(temp, src, width, height, filter, W) ;
	econvolve(dst, temp, height, width, filter, W) ;
}
	
SIFT::SIFT()
{
	O = -1;
	S = 3;
	sigman = 0.5;
	sigma0 = 1.6 * powf(2.0f, 1.0f / S) ;
	omin = -1;
	smin = -1;
	smax = S+1;
	magnif = 3.0f;
	normalizeDescriptor = true;
	temp = NULL;
	octaves = NULL;
	filter = NULL;
	pgmBuffer.data = NULL;
}
	
SIFT::SIFT(float _sigman, float _sigma0, int _O, int _S, int _omin, int _smin, int _smax)
{
	sigman = _sigman;
	sigma0 = _sigma0;
	O = _O;
	S = _S;
	omin = _omin;
	smin = _smin;
	smax= _smax;
	magnif = 3.0f;
	normalizeDescriptor = true;
	temp = NULL;
	octaves = NULL;
	filter = NULL;
	pgmBuffer.data = NULL;
}
	
	
SIFT::~SIFT()
{
	freeBuffers() ;
}

void SIFT::processImage(Img *im, std::vector<keypt_t> &keyInfo, std::vector<short> &keyDescr)
{
	printf("Processing image...\n");
	extractPgm(im, pgmBuffer);
	if(O < 1) 
	{
		O = std::max(int(std::floor(log2(std::min(pgmBuffer.width,pgmBuffer.height))) - omin -3), 1) ;
	}
	process(pgmBuffer.data, pgmBuffer.width, pgmBuffer.height);
	float threshold = 0.02;
	float edgeThreshold = 7.5;
	detectKeypoints(threshold, edgeThreshold); 
	for(KeypointsConstIter iter = keypointsBegin(); iter != keypointsEnd() ; ++iter ) 
	{
	    // detect orientations
	    float angles [4] ;
	    int nangles = computeKeypointOrientations(angles, *iter) ;

		// compute descriptors
	    for(int a = 0 ; a < nangles ; ++a) 
		{
			float descr_pt [128] ;
			computeKeypointDescriptor(descr_pt, *iter, angles[a]) ;
			for(int i=0; i<128; i++)
			{
				keyDescr.push_back((short)(descr_pt[i] * 255.0));
			}
			keyInfo.push_back(*iter);
			keyInfo[keyInfo.size()-1].orient = angles[a];
	    }
	}
	freeBuffers();
	delete [] pgmBuffer.data;
	pgmBuffer.data = NULL;
}

void SIFT::prepareImage(Img *im)
{
	printf("Processing image...\n");
	extractPgm(im, pgmBuffer);
	if(O < 1) 
	{
		O = std::max(int(std::floor(log2(std::min(pgmBuffer.width,pgmBuffer.height))) - omin -3), 1) ;
	}
	process(pgmBuffer.data, pgmBuffer.width, pgmBuffer.height);
}

void SIFT::cleanImage()
{
	freeBuffers();
	delete [] pgmBuffer.data;
	pgmBuffer.data = NULL;
}

void SIFT::processUprightImage(Img *im, std::vector<keypt_t> &keyInfo, std::vector<short> &keyDescr)
{
	printf("Processing image...\n");
	float epsilon = 1e-2;
	extractPgm(im, pgmBuffer);
	if(O < 1) 
	{
		O = std::max(int(std::floor(log2(std::min(pgmBuffer.width,pgmBuffer.height))) - omin -3), 1) ;
	}
	process(pgmBuffer.data, pgmBuffer.width, pgmBuffer.height);
	float threshold = 0.02;
	float edgeThreshold = 7.5;
	detectKeypoints(threshold, edgeThreshold); 
	for(KeypointsConstIter iter = keypointsBegin(); iter != keypointsEnd() ; ++iter ) 
	{
		// compute descriptors
		float descr_pt [128] ;
		computeKeypointDescriptor(descr_pt, *iter, 0.0) ;
		for(int i=0; i<128; i++)
		{
			keyDescr.push_back((short)(descr_pt[i] * 255.0));
		}
		keyInfo.push_back(*iter);
	}
	
	freeBuffers();
	delete [] pgmBuffer.data;
	pgmBuffer.data = NULL;
}

void SIFT::extractPgm(Img *im, PgmBuffer& buffer)
{
	buffer.width = im->width();
	buffer.height = im->height();
	if(buffer.data != NULL)
		delete [] buffer.data;
	
	buffer.data = new pixel_t[buffer.width*buffer.height];
	
	for(int y=0; y<buffer.height; y++)
	{
		for(int x=0; x<buffer.width; x++)
		{
			Color c = (*im)(x,y);
			buffer.data[y*buffer.width+x] = c[0]*0.3 + c[1]*0.59 + c[2]*0.11;
		}
	}
}
	
void SIFT::prepareBuffers()
{
	// compute buffer size
	int w = (omin >= 0) ? (width  >> omin) : (width  << -omin) ;
	int h = (omin >= 0) ? (height >> omin) : (height << -omin) ;
	int size = w*h* std::max((smax - smin), 2*((smax+1) - (smin-2) +1)) ;
		
	if( temp && tempReserved == size ) return ;
		
	freeBuffers() ;
		
	// allocate
	temp           = new pixel_t [ size ] ; 
	tempReserved   = size ;
	tempIsGrad     = false ;
	tempOctave     = 0 ;
		
	octaves = new pixel_t* [ O ] ;
	for(int o = 0 ; o < O ; ++o)
	{
		octaves[o] = new pixel_t [ (smax - smin + 1) * w * h ] ;
		w >>= 1 ;
		h >>= 1 ;
	}
}
	
void SIFT::freeBuffers()
{
	if( filter != NULL)
	{
		delete [] filter ;
	}
	filter = NULL ;
		
	if( octaves != NULL) 
	{
		for(int o = 0 ; o < O ; ++o) 
		{
			delete [] octaves[ o ] ;
		}
		delete [] octaves ;
	}
		
	octaves = NULL ;
		
	if( temp != NULL)
	{
		delete [] temp ;   
	}
	temp = NULL  ; 
}
	
keypt_t SIFT::getKeypoint(float x, float y, float sigma)
{
	/*
	The formula linking the keypoint scale sigma to the octave and
	scale index is
		 
	(1) sigma(o,s) = sigma0 2^(o+s/S)
		 
	for which
		 
	(2) o + s/S = log2 sigma/sigma0 == phi.
		 
	In addition to the scale index s (which can be fractional due to
	scale interpolation) a keypoint has an integer scale index is too
	(which is the index of the scale level where it was detected in
	the DoG scale space). We have the constraints:
		 
	- o and is are integer
		 
	- is is in the range [smin+1, smax-2  ]
		 
	- o  is in the range [omin,   omin+O-1]
		 
	- is = rand(s) most of the times (but not always, due to the way s
	is obtained by quadratic interpolation of the DoG scale space).
		
	Depending on the values of smin and smax, often (2) has multiple
	solutions is,o that satisfy all constraints.  In this case we
	choose the one with biggest index o (this saves a bit of
	computation).
		 
	DETERMINING THE OCTAVE INDEX O
		 
	From (2) we have o = phi - s/S and we want to pick the biggest
	possible index o in the feasible range. This corresponds to
	selecting the smallest possible index s. We write s = is + ds
	where in most cases |ds|<.5 (but in general |ds|<1). So we have
		
	o = phi - s/S,   s = is + ds ,   |ds| < .5 (or |ds| < 1).
		 
	Since is is in the range [smin+1,smax-2], s is in the range
	[smin+.5,smax-1.5] (or [smin,smax-1]), the number o is an integer
	in the range phi+[-smax+1.5,-smin-.5] (or
	phi+[-smax+1,-smin]). Thus the maximum value of o is obtained for
	o = floor(phi-smin-.5) (or o = floor(phi-smin)).
		 
	Finally o is clamped to make sure it is contained in the feasible
	range.
		 
	DETERMINING THE SCALE INDEXES S AND IS
		
	Given o we can derive is by writing (2) as
		 
	s = is + ds = S(phi - o).
		 
	We then take is = round(s) and clamp its value to be in the
	feasible range.
	*/
		
	int o,ix,iy,is ;
	float s,phi ;
		
	phi = log2(sigma/sigma0) ;
	o   = fast_floor( phi -  (float(smin)+.5)/S ) ;
	o   = std::min(o, omin+O-1) ;
	o   = std::max(o, omin    ) ;
	s   = S * (phi - o) ;
		
	is  = int(s + 0.5) ;
	is  = std::min(is, smax - 2) ;
	is  = std::max(is, smin + 1) ;
		
	float per = getOctaveSamplingPeriod(o) ;
	ix = int(x / per + 0.5) ;
	iy = int(y / per + 0.5) ;
		
	keypt_t key ;
	key.o  = o ;
		
	key.ix = ix ;
	key.iy = iy ;
	key.is = is ;
	
	key.x = x ;
	key.y = y ;
	key.s = s ;
		
	key.scale = sigma ;
		
	return key ;
}
	
void SIFT::process(const pixel_t* _im_pt, int _width, int _height)
{
	width  = _width ;
	height = _height ;
	prepareBuffers() ;
		
	float sigmak = powf(2.0f, 1.0 / S) ;
	float dsigma0 = sigma0 * sqrt (1.0f - 1.0f / (sigmak*sigmak) ) ;
		
	//Make pyramid base
	if( omin < 0 ) 
	{
		copyAndUpsampleRows(temp,       _im_pt, width,  height  ) ;
		copyAndUpsampleRows(octaves[0], temp,   height, 2*width ) ;      
			
		for(int o = -1 ; o > omin ; --o)
		{
			copyAndUpsampleRows(temp,       octaves[0], width  << -o,    height << -o) ;
			copyAndUpsampleRows(octaves[0], temp,       height << -o, 2*(width  << -o)) ;             
		}
			
	} 
	else if( omin > 0 ) 
	{
		copyAndDownsample(octaves[0], _im_pt, width, height, 1 << omin) ;
	} 
	else 
	{
		copy(octaves[0], _im_pt, width, height) ;
	}
		
	float sa = sigma0 * powf(sigmak, smin) ; 
	float sb = sigman / powf(2.0f,   omin) ; // review this
	if( sa > sb ) 
	{
		float sd = sqrt ( sa*sa - sb*sb ) ;
		smooth( octaves[0], temp, octaves[0], getOctaveWidth(omin), getOctaveHeight(omin), sd ) ;
	}
		
	// Make octaves
	for(int o = omin ; o < omin+O ; ++o) 
	{
		// Prepare octave base
		if( o > omin ) 
		{
			int sbest = std::min(smin + S, smax) ;
			copyAndDownsample(getLevel(o,   smin ), 
								getLevel(o-1, sbest),
								getOctaveWidth(o-1),
								getOctaveHeight(o-1), 2 ) ;
			float sa = sigma0 * powf(sigmak, smin      ) ;
			float sb = sigma0 * powf(sigmak, sbest - S ) ;
			if(sa > sb )
			{
				float sd = sqrt ( sa*sa - sb*sb ) ;
				smooth( getLevel(o,0), temp, getLevel(o,0), 
						getOctaveWidth(o), getOctaveHeight(o), sd ) ;
			}
		}
			
		// Make other levels
		for(int s = smin+1 ; s <= smax ; ++s) 
		{
			float sd = dsigma0 * powf(sigmak, s) ;
			smooth( getLevel(o,s), temp, getLevel(o,s-1),
					getOctaveWidth(o), getOctaveHeight(o), sd ) ;
		}
	}
	
}
	
void SIFT::detectKeypoints(float threshold, float edgeThreshold)
{
	keypoints.clear() ;
		
	int nValidatedKeypoints = 0 ;
		
	// Process one octave per time
	for(int o = omin ; o < omin + O ; ++o) 
	{
		int const xo = 1 ;
		int const yo = getOctaveWidth(o) ;
		int const so = getOctaveWidth(o) * getOctaveHeight(o) ;
		int const ow = getOctaveWidth(o) ;
		int const oh = getOctaveHeight(o) ;
			
		float xperiod = getOctaveSamplingPeriod(o) ;
			
		//difference of gaussians
		pixel_t* dog = temp ;
		tempIsGrad = false ;
		pixel_t* pt = dog ;
		for(int s = smin ; s <= smax-1 ; ++s) 
		{
			pixel_t* srca = getLevel(o, s  ) ;
			pixel_t* srcb = getLevel(o, s+1) ;
			pixel_t* enda = srcb ;
			while( srca != enda ) 
			{
				*pt++ = *srcb++ - *srca++ ;
			}
		}
			
		//points of extrema
		pt  = dog + xo + yo + so ;
		for(int s = smin+1 ; s <= smax-2 ; ++s)
		{
			for(int y = 1 ; y < oh - 1 ; ++y) 
			{
				for(int x = 1 ; x < ow - 1 ; ++x) 
				{          
					pixel_t v = *pt ;
					
					if( CHECK_NEIGHBORS(>,+) || CHECK_NEIGHBORS(<,-) )
					{
						keypt_t k ;
						k.ix = x ;
						k.iy = y ;
						k.is = s ;
						keypoints.push_back(k) ;
					}
					pt += 1 ;
				}
				pt += 2 ;
			}
			pt += 2*yo ;
		}
			
		// refine local max
		KeypointsIter siter ;
		KeypointsIter diter ;
				
		for(diter = siter = keypointsBegin() + nValidatedKeypoints ; siter != keypointsEnd() ; ++siter) 
		{
			int x = int( siter->ix ) ;
			int y = int( siter->iy ) ;
			int s = int( siter->is ) ;
					
			float Dx=0,Dy=0,Ds=0,Dxx=0,Dyy=0,Dss=0,Dxy=0,Dxs=0,Dys=0 ;
			float  b[3] ;
			pixel_t* pt ;
			int dx = 0 ;
			int dy = 0 ;
					
			// must be exec. at least once
			for(int iter = 0 ; iter < 5 ; ++iter)
			{
				float A[3*3] ;          
						
				x += dx ;
				y += dy ;
						
				pt = dog + xo * x + yo * y + so * (s - smin) ;
						
				#define at(dx,dy,ds) (*( pt + (dx)*xo + (dy)*yo + (ds)*so))
				// Compute the gradient
				Dx = 0.5 * (at(+1,0,0) - at(-1,0,0)) ;
				Dy = 0.5 * (at(0,+1,0) - at(0,-1,0));
				Ds = 0.5 * (at(0,0,+1) - at(0,0,-1)) ;
						
				// Compute the Hessian
				Dxx = (at(+1,0,0) + at(-1,0,0) - 2.0 * at(0,0,0)) ;
				Dyy = (at(0,+1,0) + at(0,-1,0) - 2.0 * at(0,0,0)) ;
				Dss = (at(0,0,+1) + at(0,0,-1) - 2.0 * at(0,0,0)) ;
						
				Dxy = 0.25 * ( at(+1,+1,0) + at(-1,-1,0) - at(-1,+1,0) - at(+1,-1,0) ) ;
				Dxs = 0.25 * ( at(+1,0,+1) + at(-1,0,-1) - at(-1,0,+1) - at(+1,0,-1) ) ;
				Dys = 0.25 * ( at(0,+1,+1) + at(0,-1,-1) - at(0,-1,+1) - at(0,+1,-1) ) ;
						
				// Solve linear system
				Aat(0,0) = Dxx ;
				Aat(1,1) = Dyy ;
				Aat(2,2) = Dss ;
				Aat(0,1) = Aat(1,0) = Dxy ;
				Aat(0,2) = Aat(2,0) = Dxs ;
				Aat(1,2) = Aat(2,1) = Dys ;
						
				b[0] = - Dx ;
				b[1] = - Dy ;
				b[2] = - Ds ;
						
				// Gauss elimination
				for(int j = 0 ; j < 3 ; ++j) 
				{
					// look for leading pivot
					float maxa = 0 ;
					float maxabsa = 0 ;
					int   maxi = -1 ;
					int i ;
					for(i = j ; i < 3 ; ++i) 
					{
						float a    = Aat(i,j) ;
						float absa = fabsf( a ) ;
						if ( absa > maxabsa ) 
						{
							maxa    = a ;
							maxabsa = absa ;
							maxi    = i ;
						}
					}
							
					// singular?
					if( maxabsa < 1e-10f ) 
					{
						b[0] = 0 ;
						b[1] = 0 ;
						b[2] = 0 ;
						break ;
					}
							
					i = maxi ;
							
					// swap j-th row with i-th row and
					// normalize j-th row
					for(int jj = j ; jj < 3 ; ++jj) 
					{
						std::swap( Aat(j,jj) , Aat(i,jj) ) ;
						Aat(j,jj) /= maxa ;
					}
					std::swap( b[j], b[i] ) ;
					b[j] /= maxa ;
							
					// elimination
					for(int ii = j+1 ; ii < 3 ; ++ii)
					{
						float x = Aat(ii,j) ;
						for(int jj = j ; jj < 3 ; ++jj) 
						{
							Aat(ii,jj) -= x * Aat(j,jj) ;                
						}
						b[ii] -= x * b[j] ;
					}
				}
						
				// backward substitution
				for(int i = 2 ; i > 0 ; --i) 
				{
					float x = b[i] ;
					for(int ii = i-1 ; ii >= 0 ; --ii) 
					{
						b[ii] -= x * Aat(ii,i) ;
					}
				}
						
				// If the translation of the keypoint is big, move the keypoint
				// and re-iterate the computation. Otherwise we are all set.
				dx= ((b[0] >  0.6 && x < ow-2) ?  1 : 0 )+ ((b[0] < -0.6 && x > 1   ) ? -1 : 0 ) ;
					
				dy= ((b[1] >  0.6 && y < oh-2) ?  1 : 0 )+ ((b[1] < -0.6 && y > 1   ) ? -1 : 0 ) ;
						
				if( dx == 0 && dy == 0 ) 
					break ;
			}
					
			// Accept-reject keypoint
			float val = at(0,0,0) + 0.5 * (Dx * b[0] + Dy * b[1] + Ds * b[2]) ; 
			float score = (Dxx+Dyy)*(Dxx+Dyy) / (Dxx*Dyy - Dxy*Dxy) ; 
			float xn = x + b[0] ;
			float yn = y + b[1] ;
			float sn = s + b[2] ;
						
			if(fast_abs(val) > threshold &&
			   score < (edgeThreshold+1)*(edgeThreshold+1)/edgeThreshold && 
			   score >= 0 &&
			   fast_abs(b[0]) < 1.5 &&
			   fast_abs(b[1]) < 1.5 &&
			   fast_abs(b[2]) < 1.5 &&
			   xn >= 0    &&
			   xn <= ow-1 &&
			   yn >= 0    &&
			   yn <= oh-1 &&
			   sn >= smin &&
			   sn <= smax ) 
			{
							
				diter->o  = o ;
							
				diter->ix = x ;
				diter->iy = y ;
				diter->is = s ;
							
				diter->x = xn * xperiod ; 
				diter->y = yn * xperiod ; 
				diter->s = sn ;
							
				diter->scale = getScaleFromIndex(o,sn) ;
							
				++diter ;
			}
		} // next candidate keypoint
				
		// prepare for next octave
		keypoints.resize( diter - keypoints.begin() ) ;
		nValidatedKeypoints = keypoints.size() ;
			
	} // next octave
}
	
void SIFT::prepareGrad(int o)
{ 
	int const ow = getOctaveWidth(o) ;
	int const oh = getOctaveHeight(o) ;
	int const xo = 1 ;
	int const yo = ow ;
	int const so = oh*ow ;
		
	if( ! tempIsGrad || tempOctave != o ) 
	{
		// compute dx/dy
		for(int s = smin+1 ; s <= smax-2 ; ++s) 
		{
			for(int y = 1 ; y < oh-1 ; ++y ) 
			{
				pixel_t* src  = getLevel(o, s) + xo + yo*y ;        
				pixel_t* end  = src + ow - 1 ;
				pixel_t* grad = 2 * (xo + yo*y + (s - smin -1)*so) + temp ;
				while(src != end) 
				{
					float Gx = 0.5 * ( *(src+xo) - *(src-xo) ) ;
					float Gy = 0.5 * ( *(src+yo) - *(src-yo) ) ;
					float m = fast_sqrt( Gx*Gx + Gy*Gy ) ;
					float t = fast_mod_2pi( fast_atan2(Gy, Gx) + float(2*M_PI) );
					*grad++ = pixel_t( m ) ;
					*grad++ = pixel_t( t ) ;
					++src ;
				}
			}
		}
	}
		
	tempIsGrad = true ;
	tempOctave = o ;
}
	
int	SIFT::computeKeypointOrientations(float angles [4], keypt_t keypoint)
{
	int const   nbins = 36 ;
	float const winFactor = 1.5 ;
	float hist [nbins] ;
		
	// octave
	int o = keypoint.o ;
	float xperiod = getOctaveSamplingPeriod(o) ;
		
	// offsets to move in the Gaussian scale space octave
	const int ow = getOctaveWidth(o) ;
	const int oh = getOctaveHeight(o) ;
	const int xo = 2 ;
	const int yo = xo * ow ;
	const int so = yo * oh ;
		
	// keypoint fractional geometry
	float x     = keypoint.x / xperiod ;
	float y     = keypoint.y / xperiod ;
	float sigma = keypoint.scale / xperiod ;
		
	// shall we use keypoints.ix,iy,is here?
	int xi = ((int) (x+0.5)) ; 
	int yi = ((int) (y+0.5)) ;
	int si = keypoint.is ;
		
	float_t const sigmaw = winFactor * sigma ;
	int W = (int) floor(3.0 * sigmaw) ;
		
	// skip the keypoint if it is out of bounds
	if( o  < omin   ||
		o  >=omin+O ||
		xi < 0      || 
		xi > ow-1   || 
		yi < 0      || 
		yi > oh-1   || 
		si < smin+1 || 
		si > smax-2 )
	{
		return 0 ;
	}
		
	// make sure that the gradient buffer is filled with octave o
	prepareGrad(o) ;
		
	// clear the SIFT histogram
	std::fill(hist, hist + nbins, 0) ;
		
	// fill the SIFT histogram
	pixel_t* pt = temp + xi * xo + yi * yo + (si - smin -1) * so ;
		
	#undef at
	#define at(dx,dy) (*(pt + (dx)*xo + (dy)*yo))
		
	for(int ys = std::max(-W, 1-yi) ; ys <= std::min(+W, oh -2 -yi) ; ++ys)
	{
		for(int xs = std::max(-W, 1-xi) ; xs <= std::min(+W, ow -2 -xi) ; ++xs) 
		{
			float dx = xi + xs - x;
			float dy = yi + ys - y;
			float r2 = dx*dx + dy*dy ;
				
			// limit to a circular window
			if(r2 >= W*W+0.5) 
				continue ;
				
			float wgt = fast_expn( r2 / (2*sigmaw*sigmaw) ) ;
			float mod = *(pt + xs*xo + ys*yo) ;
			float ang = *(pt + xs*xo + ys*yo + 1) ;
				
			int bin = (int) floor( nbins * ang / (2*M_PI) ) ;
			hist[bin] += mod * wgt ;        
		}
	}
		
	// smooth the histogram
	for (int iter = 0; iter < 6; iter++) 
	{
		float prev  = hist[nbins-1] ;
		float first = hist[0] ;
		int i ;
		for (i = 0; i < nbins - 1; i++)
		{
			float newh = (prev + hist[i] + hist[(i+1) % nbins]) / 3.0;
			prev = hist[i] ;
			hist[i] = newh ;
		}
		hist[i] = (prev + hist[i] + first)/3.0 ;
	}
	
	// find the histogram maximum
	float maxh = * std::max_element(hist, hist + nbins) ;
		
	// find peaks within 80% from max
	int nangles = 0 ;
	for(int i = 0 ; i < nbins ; ++i) 
	{
		float h0 = hist [i] ;
		float hm = hist [(i-1+nbins) % nbins] ;
		float hp = hist [(i+1+nbins) % nbins] ;
			
		// is this a peak?
		if( h0 > 0.8*maxh && h0 > hm && h0 > hp )
		{
			// quadratic interpolation
			float di = -0.5 * (hp - hm) / (hp+hm-2*h0) ; 
			float th = 2*M_PI * (i+di+0.5) / nbins ;      
			angles [ nangles++ ] = th ;
			if( nangles == 4 )
				goto enough_angles ;
		}
	}
	enough_angles:
		return nangles ;
}
	
//computeKeypointDescriptor()
void SIFT::normalizeHistogram(float* L_begin, float* L_end)
{
	float* L_iter ;
	float norm = 0.0 ;
			
	for(L_iter = L_begin; L_iter != L_end ; ++L_iter)
		norm += (*L_iter) * (*L_iter) ;
			
	norm = fast_sqrt(norm) ;
			
	for(L_iter = L_begin; L_iter != L_end ; ++L_iter)
		*L_iter /= (norm + std::numeric_limits<float>::epsilon() ) ;
}

void SIFT::computeKeypointDescriptor(float* descrPt, keypt_t keypoint, float angle0)
{
	//printf("scale:%d\n", keypoint.o);
	// octave
	int o = keypoint.o ;
	float xperiod = getOctaveSamplingPeriod(o) ;
		
	// offsets to move in Gaussian scale space octave
	const int ow = getOctaveWidth(o) ;
	const int oh = getOctaveHeight(o) ;
	const int xo = 2 ;
	const int yo = xo * ow ;
	const int so = yo * oh ;
		
	// keypoint fractional geometry
	float x     = keypoint.x / xperiod;
	float y     = keypoint.y / xperiod ;
	float sigma = keypoint.scale / xperiod ;
		
	float st0   = sinf( angle0 ) ;
	float ct0   = cosf( angle0 ) ;
		
	// shall we use keypoints.ix,iy,is here?
	int xi = ((int) (x+0.5)) ; 
	int yi = ((int) (y+0.5)) ;
	int si = keypoint.is ;
		
	const int NBO = 8 ;
	const int NBP = 4 ;
	const float SBP = magnif * sigma ;
	const int   W = (int) floor (sqrt(2.0) * SBP * (NBP + 1) / 2.0 + 0.5) ;
		
	// Offsets to move in the descriptor.
	// Use Lowe's convention.
	const int binto = 1 ;
	const int binyo = NBO * NBP ;
	const int binxo = NBO ;
	
	int bin ;
		
	// check bounds
	if( o  < omin   ||
		o  >=omin+O ||
		xi < 0      || 
		xi > ow-1   || 
		yi < 0      || 
		yi > oh-1   ||
		si < smin+1 ||
		si > smax-2 )
		return ;
		
	// make sure gradient buffer is up-to-date
	prepareGrad(o) ;
		
	std::fill( descrPt, descrPt + NBO*NBP*NBP, 0 ) ;
		
	//Center the scale space and the descriptor on the current keypoint. 
	//Note that dpt is pointing to the bin of center (SBP/2,SBP/2,0).

	float const * pt = temp + xi*xo + yi*yo + (si - smin - 1)*so ;
	float *dpt = descrPt + (NBP/2) * binyo + (NBP/2) * binxo ;
	
	//Process pixels in the intersection of the image rectangle
	//(1,1)-(M-1,N-1) and the keypoint bounding box.
	
	for(int dyi = std::max(-W, 1-yi) ; dyi <= std::min(+W, oh-2-yi) ; ++dyi) 
	{
		for(int dxi = std::max(-W, 1-xi) ; dxi <= std::min(+W, ow-2-xi) ; ++dxi) 
		{
			//retrieve 
			float mod   = *( pt + dxi*xo + dyi*yo + 0 ) ;
			float angle = *( pt + dxi*xo + dyi*yo + 1 ) ;
			float theta = fast_mod_2pi(-angle + angle0) ; // lowe compatible ?
				
			// fractional displacement
			float dx = xi + dxi - x;
			float dy = yi + dyi - y;
				
			// get the displacement normalized w.r.t. the keypoint
			// orientation and extension.
			float nx = ( ct0 * dx + st0 * dy) / SBP ;
			float ny = (-st0 * dx + ct0 * dy) / SBP ; 
			float nt = NBO * theta / (2*M_PI) ;
				
			// Get the gaussian weight of the sample. The gaussian window
			// has a standard deviation equal to NBP/2. Note that dx and dy
			// are in the normalized frame, so that -NBP/2 <= dx <= NBP/2.
			float const wsigma = NBP/2 ;
			float win = fast_expn((nx*nx + ny*ny)/(2.0 * wsigma * wsigma)) ;
				
			// The sample will be distributed in 8 adjacent bins.
			// We start from the ``lower-left'' bin.
			int binx = fast_floor( nx - 0.5 ) ;
			int biny = fast_floor( ny - 0.5 ) ;
			int bint = fast_floor( nt ) ;
			float rbinx = nx - (binx+0.5) ;
			float rbiny = ny - (biny+0.5) ;
			float rbint = nt - bint ;
			int dbinx ;
			int dbiny ;
			int dbint ;
				
			// Distribute the current sample into the 8 adjacent bins
			for(dbinx = 0 ; dbinx < 2 ; ++dbinx)
			{
				for(dbiny = 0 ; dbiny < 2 ; ++dbiny)
				{
					for(dbint = 0 ; dbint < 2 ; ++dbint) 
					{
						if( binx+dbinx >= -(NBP/2) &&
							binx+dbinx <   (NBP/2) &&
							biny+dbiny >= -(NBP/2) &&
							biny+dbiny <   (NBP/2) ) 
						{
							float weight = win * mod * fast_abs (1 - dbinx - rbinx) * fast_abs (1 - dbiny - rbiny) * fast_abs (1 - dbint - rbint) ;
								
							atd(binx+dbinx, biny+dbiny, (bint+dbint) % NBO) += weight ;
						}
					}            
				}
			}
		}  
	}
		
	//Standard SIFT descriptors are normalized, truncated and normalized again
	if( normalizeDescriptor ) 
	{
		//Normalize the histogram to L2 unit length.    
		normalizeHistogram(descrPt, descrPt + NBO*NBP*NBP) ;
			
		// Truncate at 0.2
		for(bin = 0; bin < NBO*NBP*NBP ; ++bin) 
		{
			if (descrPt[bin] > 0.2) 
				descrPt[bin] = 0.2;
		}
			
		// Normalize again.
		normalizeHistogram(descrPt, descrPt + NBO*NBP*NBP) ;
	}
		
}
