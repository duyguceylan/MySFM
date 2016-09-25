/*
 *  Image.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _IMG_H
#define _IMG_H

#include "Color.h"
#include <cstdio>
#include <cstring>
#include <cassert>
#include <vector>
#include <string>
#include <jpeglib.h>
#include <png.h>

using namespace std;

class Img 
{
private:
	std::vector<Color> pix;
	
public:
	int w, h;
	bool is16bit;
	
	// Constructors
	Img() : w(0), h(0), is16bit(false) {}
	Img(int w_, int h_) : pix(w_*h_), w(w_), h(h_), is16bit(false) {}
	Img(int w_, int h_, const Color &col) : pix(w_*h_, col), w(w_), h(h_), is16bit(false) {}
	Img(int w_, int h_, float val) : pix(w_*h_, Color(val)), w(w_), h(h_), is16bit(false) {}
	Img(int w_, int h_, const float *data, bool gray=false);
	Img(int w_, int h_, const std::vector<float> &data, bool gray=false);
	
	~Img() {pix.clear();}
	// Using default copy constructor, assignment operator, destructor
	
	// Array access
	float* getColors() {return &(pix[0][0]); }
	const Color &operator [] (int i) const { return pix[i]; }
	Color &operator [] (int i) { return pix[i]; }
	
	// Array access by row/column
	const Color &operator () (int x, int y) const { return pix[x + y * w]; }
	Color &operator () (int x, int y){ return pix[x + y * w]; }
	
	void setColor(int x, int y, Color c) {pix[x+y*w] = c; }
	
	// Interpolated access
	const Color lerp(float x, float y) const;
	
	// Member operators
	Img &operator += (const Img &x);
	Img &operator += (const Color &c);
	Img &operator += (const float x);
	Img &operator -= (const Img &x);
	Img &operator -= (const Color &c);
	Img &operator -= (const float x);
	Img &operator *= (const Img &x);
	Img &operator *= (const Color &c);
	Img &operator *= (const float x);
	Img &operator /= (const Img &x);
	Img &operator /= (const Color &c);
	Img &operator /= (const float x);
	Img &minimize(const Img &x);
	Img &minimize(const Color &c);
	Img &minimize(const float x);
	Img &maximize(const Img &x);
	Img &maximize(const Color &c);
	Img &maximize(const float x);
	
	// Partial compatibility with vectors
	typedef Color value_type;
	typedef Color *pointer;
	typedef const Color *const_pointer;
	typedef Color *iterator;
	typedef const Color *const_iterator;
	typedef Color &reference;
	typedef const Color &const_reference;
	typedef size_t size_type;
	typedef std::ptrdiff_t difference_type;
	
	size_t width() const {return w; }
	size_t height() const {return h; }
	size_t size() const { return w*h; }
	Color *begin(){ return &(pix[0]); }
	const Color *begin() const { return &(pix[0]); }
	Color *end() { return begin() + w*h; }
	const Color *end() const { return begin() + w*h; }
	void clear(){ pix.clear(); w=h=0; is16bit = false; }
	bool empty() const { return (w <= 0 || h <= 0); }
	void resize(int w_, int h_) { w = w_; h = h_; pix.resize(w*h); }
	void resize(int w_, int h_, const Color &col) { w = w_; h = h_; pix.resize(w*h, col); }
	void resize(int w_, int h_, float val) { w = w_; h = h_; pix.resize(w*h, Color(val,val,val)); }
	const Color sum() const;
	const Color avg() const { return sum() / float(w*h); }
	const Color minColor() const;
	const Color maxColor() const;
	
	// Color transformations.  These are done in-place, unlike
	// the functions on Colors.
	void col_transform(float m11, float m12, float m13,
					   float m21, float m22, float m23,
					   float m31, float m32, float m33);
	void convert(Color::Colorspace src, Color::Colorspace dst);
	void gamma(float g), gamma(Color::Colorspace dst);
	void ungamma(float g), ungamma(Color::Colorspace dst);
	
	// Simple transformations
	void flipX(), flipY(), rotCW(), rotCCW();
	void crop(int start_x, int start_y, int new_w, int new_h);
	Img* getCroppedImg(int start_x, int start_y, int new_w, int new_h);
	void set(int start_x, int start_y, int subimg_w, int subimg_h, const Color &c);
	void set(const Color &c) { set(0, 0, w, h, c); }
	void set(int start_x, int start_y, int subimg_w, int subimg_h, const Img &im2);
	void set(int start_x, int start_y, const Img &im2) { set(start_x, start_y, im2.w, im2.h, im2); }
	
	void computeIntegralImage(vector<float> &intergralImg, vector<float> &sqIntegralImg);
	
	Color getDivergenceAtPixel(int x, int y);
	
	// Input/output
	bool read(const std::string &filename);
	bool write(const std::string &filename);
	
	void gaussianSmooth(float sigma, float **smoothedImgR, float **smoothedImgG, float **smoothedImgB);
	void downScaleAndGaussianSmoothImage();
	float computeNCCScoreWithGrayScale(Img neighborImg);
private:
	bool read_pbm(std::FILE *f);
	bool read_jpg(std::FILE *f);
	bool read_png(std::FILE *f);
	bool write_jpg(std::FILE *f);
	bool write_png(std::FILE *f);
	static inline bool we_are_little_endian()
	{ 
		int tmp = 1;
		return !!(* (unsigned char *) &tmp); 
	}
};

// Nonmember operators
static inline const Img operator + (const Img &x, const Img &y) { return Img(x) += y; }
static inline const Img operator - (const Img &x, const Img &y) { return Img(x) -= y; }
static inline const Img operator * (const Img &x, const Img &y) { return Img(x) *= y; }
static inline const Img operator * (const Img &x, const float y) { return Img(x) *= y; }
static inline const Img operator * (const float x, const Img &y) { return y * x; }
static inline const Img operator * (const Img &x, const Color &c) { return Img(x) *= c; }
static inline const Img operator * (const Color &c, const Img &y) { return y * c; }
static inline const Img operator / (const Img &x, const Img &y) { return Img(x) /= y; }
static inline const Img operator / (const Img &x, const float y) { return Img(x) /= y; }
static inline const Img operator / (const Img &x, const Color &c) { return Img(x) /= c; }
static inline const Img operator / (const Color &c, const Img &y)
{
	Img result(y);
	int n = result.w * result.h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		result[i] = c / result[i];
	return result;
}

static inline const Img &operator + (const Img &x) { return x; }
static inline const Img operator - (const Img &x)
{
	Img result(x);
	int n = result.w * result.h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		result[i] = -result[i];
	return result;
}

static inline bool operator ! (const Img &x) { return x.empty(); }

// Other nonmember functions
static inline const Img minimizeImage(const Img &x, const Img &y) { return Img(x).minimize(y); }
static inline const Img maximizeImage(const Img &x, const Img &y) { return Img(x).maximize(y); }

#endif
