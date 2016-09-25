/*
 *  Image.cpp
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "Image.h"
#include "MyMatrix.h"

float clamp(const float &x, const float &a, const float &b)
{
	return x > a ? x < b ? x : b : a;  // returns a on NaN
};

Img::Img(int w_, int h_, const float *data, bool gray)
{
	pix.resize(w_*h_);
	w = w_;
	h = h_;
	is16bit = false;
	if (gray)
		for (int i=0; i < w*h; i++) 
			pix[i] = Color(data[i]);
	else
		for (int i=0; i < w*h; i++) 
			pix[i] = Color(&data[3*i]);
}

Img::Img(int w_, int h_, const std::vector<float> &data, bool gray)
{
	pix.resize(w_*h_);
	w = w_;
	h = h_;
	is16bit = false;
	if (gray)
		for (int i=0; i < w*h; i++) 
			pix[i] = Color(data[i]);
	else
		for (int i=0; i < w*h; i++)
			pix[i] = Color(&data[3*i]);
}

const Color Img::lerp(float x, float y) const
{
	x = clamp(x, 0.0f, w - 1.0f);
	y = clamp(y, 0.0f, h - 1.0f);
	int X = int(x), Y = int(y);
	float fx = x - float(X);
	float fy = y - float(Y);
	const Color &ll = pix[X + Y * w];
	const Color &lr = (fx > 0.0f) ? pix[X+1 + Y * w] : ll;
	const Color &ul = (fy > 0.0f) ? pix[X + (Y+1) * w] : ll;
	const Color &ur = (fx > 0.0f) ? ((fy > 0.0f) ? pix[(X+1) + (Y+1) * w] : lr) : ((fy > 0.0f) ? ul : ll);
	return (1.0f - fx)*(ll*(1.0f - fy) + ul*fy) + (fx  * (lr*(1.0f - fy)) + ur*fy);
	
}

Img& Img::operator += (const Img &x)
{
	assert(w == x.w);
	assert(h == x.h);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] += x[i];
	return *this;
}

Img& Img::operator += (const Color &c)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] += c;
	return *this;
}

Img& Img::operator += (const float x)
{
	Color c(x);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] += c;
	return *this;
}

Img& Img::operator -= (const Img &x)
{
	assert(w == x.w);
	assert(h == x.h);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] -= x[i];
	return *this;
}

Img& Img::operator -= (const Color &c)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] -= c;
	return *this;
}

Img& Img::operator -= (const float x)
{
	Color c(x);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] -= c;
	return *this;
}

Img& Img::operator *= (const Img &x)
{
	assert(w == x.w);
	assert(h == x.h);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] *= x[i];
	return *this;
}

Img& Img::operator *= (const Color &c)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] *= c;
	return *this;
}

Img& Img::operator *= (const float x)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] *= x;
	return *this;
}

Img& Img::operator /= (const Img &x)
{
	assert(w == x.w);
	assert(h == x.h);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] /= x[i];
	return *this;
}

Img& Img::operator /= (const Color &c)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] /= c;
	return *this;
}

Img& Img::operator /= (const float x)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] /= x;
	return *this;
}

Img& Img::minimize(const Img &x)
{
	assert(w == x.w);
	assert(h == x.h);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i].min(x[i]);
	return *this;
}

Img& Img::minimize(const Color &c)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i].min(c);
	return *this;
}

Img& Img::minimize(const float x)
{
	Color c(x);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i].min(c);
	return *this;
}

Img& Img::maximize(const Img &x)
{
	assert(w == x.w);
	assert(h == x.h);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i].min(x[i]);
	return *this;
}

Img& Img::maximize(const Color &c)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i].max(c);
	return *this;
}

Img& Img::maximize(const float x)
{
	Color c(x);
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i].max(c);
	return *this;
}

const Color Img::sum() const 
{ 
	Color total; 
	int n = w*h;
	for (int i = 0; i < n; i++) 
		total += pix[i];
	return total; 
}

const Color Img::minColor() const
{ 
	int n = w*h;
	if (!n) 
		return Color();
	Color m = pix[0];
	for (int i = 1; i < n; i++) 
		m.min(pix[i]);
	return m; 
}

const Color Img::maxColor() const
{ 
	int n = w*h;
	if (!n)
		return Color();
	Color m = pix[0];
	for (int i = 1; i < n; i++) 
		m.max(pix[i]);
	return m; 
}

// Color transformations
void Img::col_transform(float m11, float m12, float m13,
						float m21, float m22, float m23,
						float m31, float m32, float m33)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] = pix[i].col_transform(m11,m12,m13,m21,m22,m23,m31,m32,m33);
}

inline void Img::convert(Color::Colorspace src, Color::Colorspace dst)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] = pix[i].convert(src,dst);
}

void Img::gamma(float g)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] = pix[i].gamma(g);
}

void Img::gamma(Color::Colorspace dst)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] = pix[i].gamma(dst);
}

void Img::ungamma(float g)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] = pix[i].ungamma(g);
}

void Img::ungamma(Color::Colorspace dst)
{
	int n = w*h;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		pix[i] = pix[i].ungamma(dst);
}

// Geometric transformations
void Img::flipX()
{
	int w2 = w/2;
#pragma omp parallel for
	for (int y = 0; y < h; y++) 
	{
		int row = y*w, rowend = row + w-1;
		for (int x = 0; x < w2; x++)
			std::swap(pix[row+x], pix[rowend-x]);
	}
}

void Img::flipY()
{
	int h2 = h/2;
#pragma omp parallel for
	for (int y = 0; y < h2; y++) 
	{
		int row = y*w, other = (h-1-y)*w;
		for (int x = 0; x < w; x++)
			std::swap(pix[x+row], pix[x+other]);
	}
}

void Img::rotCW()
{
	std::vector<Color> tmp(pix);
	std::swap(w,h);
#pragma omp parallel for
	for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
			pix[x+y*w] = tmp[y+(w-1-x)*h];
}

void Img::rotCCW()
{
	std::vector<Color> tmp(pix);
	std::swap(w,h);
#pragma omp parallel for
	for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
			pix[x+y*w] = tmp[h-1-y+x*h];
}

void Img::crop(int start_x, int start_y, int new_w, int new_h)
{
	std::vector<Color> tmp(pix);
	pix.resize(new_w * new_h);
#pragma omp parallel for
	for (int y = 0; y < new_h; y++)
		for (int x = 0; x < new_w; x++)
			pix[x+y*new_w] = tmp[x+start_x+(y+start_y)*w];
	w = new_w;
	h = new_h;
}

Img* Img::getCroppedImg(int start_x, int start_y, int new_w, int new_h)
{
	Img *newImg = new Img(new_w,new_h);
	#pragma omp parallel for
	for (int y = 0; y < new_h; y++)
		for (int x = 0; x < new_w; x++)
			newImg->setColor(x, y, (*this)(x+start_x,(y+start_y)));
	return newImg;
}

// Set a subimage
void Img::set(int start_x, int start_y,int subimg_w, int subimg_h, const Color &c)
{
	if (start_x >= w || start_y >= h)
		return;
	if (subimg_w <= 0 || subimg_h <= 0)
		return;
	int end_x = start_x + subimg_w;
	int end_y = start_y + subimg_h;
	if (end_x < 0 || end_y < 0)
		return;
	if (start_x < 0) start_x = 0;
	if (start_y < 0) start_y = 0;
	if (end_x > w) end_x = w;
	if (end_y > h) end_y = h;
	subimg_w = end_x - start_x;
	subimg_h = end_y - start_y;
	
#pragma omp parallel for
	for (int y = 0; y < subimg_h; y++)
		for (int x = 0; x < subimg_w; x++)
			pix[start_x+x+w*(start_y+y)] = c;
}

void Img::set(int start_x, int start_y, int subimg_w, int subimg_h, const Img &im2)
{
	if (start_x >= w || start_y >= h)
		return;
	if (subimg_w <= 0 || subimg_h <= 0)
		return;
	if (subimg_w > im2.w)
		subimg_w = im2.w;
	if (subimg_h > im2.h)
		subimg_h = im2.h;
	int end_x = start_x + subimg_w;
	int end_y = start_y + subimg_h;
	if (end_x < 0 || end_y < 0)
		return;
	if (start_x < 0) start_x = 0;
	if (start_y < 0) start_y = 0;
	if (end_x > w) end_x = w;
	if (end_y > h) end_y = h;
	subimg_w = end_x - start_x;
	subimg_h = end_y - start_y;
	
#pragma omp parallel for
	for (int y = 0; y < subimg_h; y++)
		for (int x = 0; x < subimg_w; x++)
			pix[start_x+x+w*(start_y+y)] = im2(x,y);
}

void Img::computeIntegralImage(vector<float> &integral, vector<float> &sqIntegral)
{
	int minCoord = w;
	if(h < minCoord)
		minCoord = h;
	
	for(int i=0; i<minCoord; i++)
	{
		//fill row
		for(int x=i; x<w; x++)
		{
			float v = ((*this)(x,i)[0]+(*this)(x,i)[1]+(*this)(x,i)[2]) / 3.0;
			float vSq = v*v;
			if(x > 0)
			{
				v += integral[x-1 + i*w];
				vSq += sqIntegral[x-1 + i*w];
			}
			if(i > 0)
			{
				v += integral[x+ (i-1)*w];
				vSq += sqIntegral[x + (i-1)*w];
			}
			if(x>0 && i>0)
			{
				v -= integral[x-1 + (i-1)*w];
				vSq -= sqIntegral[x-1 + (i-1)*w]; 
			}
			integral[x + i*w] = v;
			sqIntegral[x + i*w] = vSq;
		}
		//fill column
		for(int y=i; y<h; y++)
		{
			float v = ((*this)(i,y)[0]+(*this)(i,y)[1]+(*this)(i,y)[2]) / 3.0;
			float vSq = v*v;
			if(y > 0)
			{
				v += integral[i + (y-1)*w];
				vSq += sqIntegral[i + (y-1)*w];
			}
			if(i > 0)
			{
				v += integral[i-1 + y*w];
				vSq += sqIntegral[i-1 + y*w];
			}
			if(y>0 && i>0)
			{
				v -= integral[i-1 + (y-1)*w];
				vSq -= sqIntegral[i-1 + (y-1)*w];
			}
			integral[i + y*w] = v;
			sqIntegral[i + y*w] = vSq;
		}
	}
	
	if(minCoord == h)
	{
		//fill remaining columns
		for(int i=minCoord; i<w; i++)
		{
			for(int y=i; y<h; y++)
			{
				float v = ((*this)(i,y)[0]+(*this)(i,y)[1]+(*this)(i,y)[2]) / 3.0;
				float vSq = v*v;
				if(y > 0)
				{
					v += integral[i + (y-1)*w];
					vSq += sqIntegral[i + (y-1)*w];
				}
				if(i > 0)
				{
					v += integral[i-1 + y*w];
					vSq += sqIntegral[i-1 + y*w];
				}
				if(y>0 && i>0)
				{
					v -= integral[i-1 + (y-1)*w];
					vSq -= sqIntegral[i-1 + (y-1)*w];
				}
				integral[i + y*w] = v;
				sqIntegral[i + y*w] = vSq;
			}			
		}
	}
	else if(minCoord == w)
	{
		//fill remaining rows
		for(int i=minCoord; i<h; i++)
		{
			for(int x=i; x<w; x++)
			{
				float v = ((*this)(x,i)[0]+(*this)(x,i)[1]+(*this)(x,i)[2]) / 3.0;
				float vSq = v*v;
				if(x > 0)
				{
					v += integral[x-1 + i*w];
					vSq += sqIntegral[x-1 + i*w];
				}
				if(i > 0)
				{
					v += integral[x+ (i-1)*w];
					vSq += sqIntegral[x + (i-1)*w];
				}
				if(x>0 && i>0)
				{
					v -= integral[x-1 + (i-1)*w];
					vSq -= sqIntegral[x-1 + (i-1)*w]; 
				}
				integral[x + i*w] = v;
				sqIntegral[x + i*w] = vSq;
			}
		}
	}
}

// I/O
bool Img::read(const std::string &filename)
{
	using namespace std;
	
	FILE *f = strcmp(filename.c_str(), "-") ?
	fopen(filename.c_str(), "rb") : stdin;
	if (!f) 
	{
		fprintf(stderr, "Couldn't open %s\n", filename.c_str());
		return false;
	}
	
	clear();
	
	// Parse magic number
	int m1 = fgetc(f), m2 = fgetc(f);
	ungetc(m2,f); ungetc(m1,f);
	
	if (m1 == 'P' && m2 == '4')
		return read_pbm(f);
	if (m1 == 0xff && m2 == 0xd8)
		return read_jpg(f);
	if (m1 == 0x89 && m2 == 'P')
		return read_png(f);
	
	bool is_pfm = false;
	int channels = 3;
	if (m1 == 'P' && m2 == '5')
		channels = 1;
	else if (m1 == 'P' && m2 == 'F')
		is_pfm = true;
	else if (m1 == 'P' && m2 == 'f')
		is_pfm = true, channels = 1;
	else if (!(m1 == 'P' && m2 == '6')) 
	{
		fclose(f);
		fprintf(stderr, "Unknown file type\n");
		return false;
	}
	
	char buf[1024];
	fgets(buf, sizeof(buf), f);
	fgets(buf, sizeof(buf), f);
	while (buf[0] == '#')
		fgets(buf, sizeof(buf), f);
	
	// Get size
	if (sscanf(buf, "%d %d", &w, &h) != 2) 
	{
		fclose(f);
		w = h = 0;
		fprintf(stderr, "Couldn't read dimensions\n");
		return false;
	}
	fgets(buf, sizeof(buf), f);
	while (buf[0] == '#')
		fgets(buf, sizeof(buf), f);
	
	// Get maxval
	bool need_swap = Img::we_are_little_endian();
	float maxval;
	if (is_pfm) 
	{
		if (sscanf(buf, "%f", &maxval) != 1) 
		{
			fclose(f);
			w = h = 0;
			fprintf(stderr, "Couldn't read maxval\n");
			return false;
		}
		if (maxval < 0.0f) 
		{
			maxval = -maxval;
			need_swap = !need_swap;
		}
	} 
	else 
	{
		int m;
		if (sscanf(buf, "%d", &m) != 1 || m < 1 || m > 65535) 
		{
			fclose(f);
			w = h = 0;
			fprintf(stderr, "Couldn't read maxval\n");
			return false;
		}
		maxval = m;
		is16bit = (m >= 256);
	}
	
	float scale = 1.0f / maxval;
	
	// Read data
	int n = w * h;
	int nbytes = is_pfm ? 4*channels*n : (1+int(is16bit))*channels*n;
	std::vector<unsigned char> data(nbytes);
	if (!fread(&data[0], nbytes, 1, f)) 
	{
		fclose(f);
		w = h = 0;
		fprintf(stderr, "Couldn't read image pixels\n");
		return false;
	}
	
	fclose(f);
	pix.resize(n);
	
	if (is_pfm) 
	{
		if (need_swap) 
		{
#pragma omp parallel for
			for (int i = 0; i < nbytes; i += 4) 
			{
				std::swap(data[i  ], data[i+3]);
				std::swap(data[i+1], data[i+2]);
			}
		}
		const float *fdata = (const float *) &data[0];
		if (channels == 1) 
		{
#pragma omp parallel for
			for (int i = 0; i < n; i++)
				pix[i] = fdata[i];
		} 
		else 
		{
#pragma omp parallel for
			for (int i = 0; i < n; i++)
				pix[i] = Color(&fdata[3*i]);
		}
	} 
	else if (is16bit) 
	{
		if (channels == 1) 
		{
#pragma omp parallel for
			for (int i = 0; i < n; i++) 
			{
				int p = 256*data[2*i] + data[2*i+1];
				pix[i] = scale * p;
			}
		} 
		else 
		{
#pragma omp parallel for
			for (int i = 0; i < n; i++) 
			{
				int p1 = 256*data[6*i  ] + data[6*i+1];
				pix[i][0] = scale * p1;
				int p2 = 256*data[6*i+2] + data[6*i+3];
				pix[i][1] = scale * p2;
				int p3 = 256*data[6*i+4] + data[6*i+5];
				pix[i][2] = scale * p3;
			}
		}
	} 
	else 
	{ 
		// 8-bit
		if (channels == 1) 
		{
#pragma omp parallel for
			for (int i = 0; i < n; i++)
				pix[i] = scale * data[i];
		} 
		else 
		{
#pragma omp parallel for
			for (int i = 0; i < n; i++) 
			{
				pix[i][0] = scale * data[3*i  ];
				pix[i][1] = scale * data[3*i+1];
				pix[i][2] = scale * data[3*i+2];
			}
		}
	}
	
	return true;
}

bool Img::read_pbm(std::FILE *f)
{
	using namespace std;
	
	char buf[1024];
	fgets(buf, sizeof(buf), f);
	fgets(buf, sizeof(buf), f);
	while (buf[0] == '#')
		fgets(buf, sizeof(buf), f);
	
	// Get size
	if (sscanf(buf, "%d %d", &w, &h) != 2) 
	{
		fclose(f);
		w = h = 0;
		fprintf(stderr, "Couldn't read dimensions\n");
		return false;
	}
	
	// Read data
	pix.resize(w * h);
	int bytes_per_row = (w + 7) / 8;
	std::vector<unsigned char> data(bytes_per_row);
	int ind = 0;
	for (int i = 0; i < h; i++) 
	{
		if (!fread(&data[0], bytes_per_row, 1, f)) 
		{
			fclose(f);
			w = h = 0;
			pix.clear();
			fprintf(stderr, "Couldn't read image pixels\n");
			return false;
		}
		for (unsigned j = 0; j < w; j++) 
		{
			unsigned byte = j >> 3u;
			unsigned bit = 7u - (j & 7u);
			if ((data[byte] >> bit) & 1u)
				pix[ind++] = Color::black();
			else
				pix[ind++] = Color::white();
		}
	}
	
	fclose(f);
	return true;
}

bool Img::read_jpg(std::FILE *f)
{
	using namespace std;
	
	jpeg_decompress_struct cinfo;
	jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, f);
	jpeg_read_header(&cinfo, TRUE);
	cinfo.out_color_space = JCS_RGB;
	jpeg_start_decompress(&cinfo);
	w = cinfo.output_width;
	h = cinfo.output_height;
	pix.resize(w*h);
	std::vector<unsigned char> buf(3*w);
	for (int i = 0; i < h; i++) 
	{
		JSAMPROW rowptr = (JSAMPROW) &buf[0];
		jpeg_read_scanlines(&cinfo, &rowptr, 1);
		for (int j = 0; j < w; j++)
			pix[w*i+j] = Color(buf[3*j], buf[3*j+1], buf[3*j+2]);
	}
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	fclose(f);
	return true;
}

bool Img::read_png(std::FILE *f)
{
	using namespace std;
	
	png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	png_infop info_ptr = png_create_info_struct(png_ptr);
	png_infop end_ptr = png_create_info_struct(png_ptr);
	png_init_io(png_ptr, f);
	png_read_info(png_ptr, info_ptr);
	png_set_expand(png_ptr);
	png_set_strip_alpha(png_ptr);
	png_set_gray_to_rgb(png_ptr);
    
    png_uint_32 width, height;
    int bit_depth, color_type, interlace_type;
    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth,&color_type, &interlace_type, NULL, NULL);
    
	w = width;
	h = height;
    
	int n = w*h;
	pix.resize(n);
	if (bit_depth == 16) 
	{
		is16bit = true;
		if (Img::we_are_little_endian())
			png_set_swap(png_ptr);
		std::vector<unsigned short> buf(3*w*h);
		std::vector<png_bytep> row_pointers(h);
		for (int i = 0; i < h; i++)
			row_pointers[i] = (png_bytep) &buf[3*w*i];
		png_read_image(png_ptr, &row_pointers[0]);
		float scale = 1.0f / 65535;
		for (int i = 0; i < n; i++) 
		{
			pix[i][0] = scale * buf[3*i];
			pix[i][1] = scale * buf[3*i+1];
			pix[i][2] = scale * buf[3*i+2];
		}
	} 
	else 
	{
		std::vector<unsigned char> buf(3*w*h);
		std::vector<png_bytep> row_pointers(h);
		for (int i = 0; i < h; i++)
			row_pointers[i] = (png_bytep) &buf[3*w*i];
		png_read_image(png_ptr, &row_pointers[0]);
		float scale = 1.0f / 255;
		for (int i = 0; i < n; i++) 
		{
			pix[i][0] = scale * buf[3*i];
			pix[i][1] = scale * buf[3*i+1];
			pix[i][2] = scale * buf[3*i+2];
		}
	}
	png_read_end(png_ptr, end_ptr);
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_ptr);
	fclose(f);
	return true;
}

bool Img::write(const std::string &filename)
{
	using namespace std;
	
	FILE *f = strcmp(filename.c_str(), "-") ? fopen(filename.c_str(), "wb") : stdout;
	if (!f) 
	{
		fprintf(stderr, "Couldn't open %s\n", filename.c_str());
		return false;
	}
	
	const char *dot = strrchr(filename.c_str(), '.');
	if (dot && !strcasecmp(dot, ".jpg"))
		return write_jpg(f);
	else if (dot && !strcasecmp(dot, ".png"))
		return write_png(f);
	else if (dot && !strcasecmp(dot, ".pfm")) 
	{
		fprintf(f, "PF\n%d %d\n%.1f\n", w, h,Img::we_are_little_endian() ? -1.0f : 1.0f);
		fwrite(&pix[0][0], 4*3*w*h, 1, f);
		fclose(f);
		return true;
	}
	
	// else write PPM
	fprintf(f, "P6\n%d %d\n%d\n", w, h, is16bit ? 65535 : 255);
	if (is16bit) 
	{
		for (int i = 0; i < w*h; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				float p = pix[i][j];
				unsigned short s = (unsigned short) clamp(int(round(p * 65535.0f)), 0, 65535);
				unsigned char c[2] = { (s >> 8u), (s & 0xffu) };
				fwrite(&c, 2, 1, f);
			}
		}
	} 
	else 
	{
		for (int i = 0; i < w*h; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				float p = pix[i][j];
				unsigned char c = (unsigned char) clamp(int(round(p * 255.0f)), 0, 255);
				fwrite(&c, 1, 1, f);
			}
		}
	}
	fclose(f);
	return true;
}

bool Img::write_jpg(std::FILE *f)
{
	using namespace std;
	
	jpeg_compress_struct cinfo;
	jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	cinfo.image_width = w;
	cinfo.image_height = h;
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;
	jpeg_set_defaults(&cinfo);  
	jpeg_set_quality(&cinfo, 90, TRUE);
	cinfo.optimize_coding = TRUE;
	cinfo.dct_method = JDCT_FLOAT;
	cinfo.comp_info[0].h_samp_factor = 1;
	cinfo.comp_info[0].v_samp_factor = 1;
	jpeg_stdio_dest(&cinfo, f);
	jpeg_start_compress(&cinfo, TRUE);
	std::vector<unsigned char> buf(3*w);
	JSAMPROW rowptr = (JSAMPROW) &buf[0];
	for (int i = 0; i < h; i++) 
	{
		for (int j = 0; j < w; j++) 
		{
			const Color &p = pix[w*i+j];
			buf[3*j  ] = clamp(int(round(p[0] * 255.0f)), 0, 255);
			buf[3*j+1] = clamp(int(round(p[1] * 255.0f)), 0, 255);
			buf[3*j+2] = clamp(int(round(p[2] * 255.0f)), 0, 255);
		}
		jpeg_write_scanlines(&cinfo, &rowptr, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	fclose(f);
	return true;
}

bool Img::write_png(std::FILE *f)
{
	using namespace std;
	
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	png_infop info_ptr = png_create_info_struct(png_ptr);
	png_init_io(png_ptr, f);
	png_set_IHDR(png_ptr, info_ptr, w, h, is16bit ? 16 : 8, 
				 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
				 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);
	if (is16bit) 
	{
		if (Img::we_are_little_endian())
			png_set_swap(png_ptr);
		std::vector<unsigned short> buf(3*w);
		for (int i = 0; i < h; i++) 
		{
			for (int j = 0; j < w; j++) 
			{
				const Color &p = pix[w*i+j];
				buf[3*j  ] = clamp(int(round(p[0] * 65535.0f)), 0, 65535);
				buf[3*j+1] = clamp(int(round(p[1] * 65535.0f)), 0, 65535);
				buf[3*j+2] = clamp(int(round(p[2] * 65535.0f)), 0, 65535);
			}
			png_write_row(png_ptr, (png_bytep) &buf[0]);
		}
	} 
	else 
	{
		std::vector<unsigned char> buf(3*w);
		for (int i = 0; i < h; i++) 
		{
			for (int j = 0; j < w; j++) 
			{
				const Color &p = pix[w*i+j];
				buf[3*j  ] = clamp(int(round(p[0] * 255.0f)), 0, 255);
				buf[3*j+1] = clamp(int(round(p[1] * 255.0f)), 0, 255);
				buf[3*j+2] = clamp(int(round(p[2] * 255.0f)), 0, 255);
			}
			png_write_row(png_ptr, (png_bytep) &buf[0]);
		}
	}
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	fclose(f);
	return true;
}

void makeGaussianKernel(float sigma, float **kernel, int *windowsize)
{
	int i, center;
	float x, fx, sum=0.0;
	
	*windowsize = 1 + 2 * ceil(2.5 * sigma);
	center = (*windowsize) / 2;
	
	(*kernel) = new float[*windowsize];
	
	for(i=0;i<(*windowsize);i++)
	{
		x = (float)(i - center);
		fx = exp(-0.5*x*x/(sigma*sigma)) / (sigma * sqrt(2.0*M_PI));
		(*kernel)[i] = fx;
		sum += fx;
	}
	
	for(i=0;i<(*windowsize);i++) 
		(*kernel)[i] /= sum;
}

void Img::gaussianSmooth(float sigma, float **smoothedImgR, float **smoothedImgG, float **smoothedImgB)
{
	int x, y, incr;
	int windowsize;
	float *tempImgR;
	float *tempImgG;
	float *tempImgB;
	float *kernel;
	float dotR, dotG, dotB, sum;
	Color color;
	
	makeGaussianKernel(sigma, &kernel, &windowsize);
	
	int center = windowsize / 2;
	
	tempImgR = new float[w*h];
	tempImgG = new float[w*h];
	tempImgB = new float[w*h];
	(*smoothedImgR) = new float[w*h];
	(*smoothedImgG) = new float[w*h];
	(*smoothedImgB) = new float[w*h];
	
	//blur in the x direction
	for(x=0;x<w;x++)
	{
		for(y=0;y<h;y++)
		{
			dotR = dotG = dotB = 0.0;
			sum = 0.0;
			for(incr=(-center);incr<=center;incr++)
			{
				if(((x+incr) >= 0) && ((x+incr) < w))
				{
					color = pix[y*w + x+incr];//[0] + pix[y*w + x+incr][1] + pix[y*w + x+incr][2]) / 3.0;
					dotR += color[0] * kernel[center+incr];
					dotG += color[1] * kernel[center+incr];
					dotB += color[2] * kernel[center+incr];
					sum += kernel[center+incr];
				}
			}
			tempImgR[x+y*w] = dotR/sum;
			tempImgG[x+y*w] = dotG/sum;
			tempImgB[x+y*w] = dotB/sum;
		}
	}
	
	//blur in the y direction
	for(y=0;y<h;y++)
	{
		for(x=0;x<w;x++)
		{
			sum = 0.0;
			dotR = dotG = dotB = 0.0;
			for(incr=(-center);incr<=center;incr++)
			{
				if(((y+incr) >= 0) && ((y+incr) < h))
				{
					dotR += tempImgR[(y+incr)*w+x] * kernel[center+incr];
					dotG += tempImgG[(y+incr)*w+x] * kernel[center+incr];
					dotB += tempImgB[(y+incr)*w+x] * kernel[center+incr];
					sum += kernel[center+incr];
				}
			}
			(*smoothedImgR)[y*w+x] = dotR/sum;
			(*smoothedImgG)[y*w+x] = dotG/sum;
			(*smoothedImgB)[y*w+x] = dotB/sum;
		}
	}
	
	delete tempImgR;
	delete tempImgG;
	delete tempImgB;
	delete kernel;
}

void Img::downScaleAndGaussianSmoothImage()
{
	//apply gaussian filter
	Matrix4f mask;
	mask[0] = Vec4f(1.0, 3.0, 3.0, 1.0);  mask[1] = Vec4f(3.0, 9.0, 9.0, 3.0);
	mask[2] = Vec4f(3.0, 9.0, 9.0, 3.0);  mask[3] = Vec4f(1.0, 3.0, 3.0, 1.0);
	
	float total = 64.0f;
	mask /= total;
	
	int width = w / 2;
	int height = h / 2;
	Img result(width, height);
		
	for (int y = 0; y < height; ++y) 
	{      
		for (int x = 0; x < width; ++x) 
		{
			Color color(0.0, 0.0, 0.0);
			Color tmpColor;
			float denom = 0.0;
				
			for (int j = -1; j < 3; ++j) 
			{
				int ytmp = 2 * y + j;
				if (ytmp < 0 || h - 1 < ytmp)
					continue;
					
				for (int i = -1; i < 3; ++i) 
				{
					int xtmp = 2 * x + i;
					if (xtmp < 0 || w - 1 < xtmp)
						continue;
						
					tmpColor = (*this)(xtmp, ytmp);
					tmpColor[0] = mask[i+1][j+1] * tmpColor[0];
					tmpColor[1] = mask[i+1][j+1] * tmpColor[1];
					tmpColor[2] = mask[i+1][j+1] * tmpColor[2];
						
					color += tmpColor;
						
					denom += mask[i+1][j+1];
				}
			}
				
			color /= denom;
			result.setColor(x,y,color);
		}
	}
	
	pix.clear();
	for (int y = 0; y < height; ++y) 
	{      
		for (int x = 0; x < width; ++x) 
		{
			pix[y*width+x] = result(x,y);
		}
	}
	w = width;
	h = height;
}

float Img::computeNCCScoreWithGrayScale(Img neighborImg)
{
	float score = 0.0;
	
	float mean1 = 0.0;
	float mean2 = 0.0;
	float stdDev1 = 0.0;
	float stdDev2 = 0.0;
	float r, n;
	
	int pixelCount = size();
	int pixelCount2 = neighborImg.size();
	if(pixelCount != pixelCount2)
	{
		printf("The number of pixels in two images don't match.\n");
		return -1.0;
	}
	
	for(int i=0; i<pixelCount; i++)
	{
		r = (pix[i][0] + pix[i][1] + pix[i][2]) / 3.0;
		mean1 += r;
		stdDev1 += r * r;
		
		n = (neighborImg[i][0] + neighborImg[i][1] + neighborImg[i][2]) / 3.0;
		mean2 += n;
		stdDev2 += n * n;
	}
	
	mean1 /= pixelCount;
	stdDev1 = stdDev1/pixelCount - mean1*mean1;
	if(stdDev1 < 0.0) stdDev1 = 0.0;
	stdDev1 = sqrt(stdDev1);
	
	
	mean2 /= pixelCount;
	stdDev2 = stdDev2/pixelCount - mean2*mean2;
	if(stdDev2 < 0.0) stdDev2 = 0.0;
	stdDev2 = sqrt(stdDev2);
	
	
	if(stdDev1 == 0.0)
		stdDev1 = 1.0;
	
	if(stdDev2 == 0.0)
		stdDev2 = 1.0;
	
	for(int i=0; i<pixelCount; i++)
	{
		r = (pix[i][0] + pix[i][1] + pix[i][2]) / 3.0;
		n = (neighborImg[i][0] + neighborImg[i][1] + neighborImg[i][2]) / 3.0;
		
		score += (r-mean1)*(n-mean2);
		if(i==(int)(pixelCount/2) && score < (pixelCount/2) * stdDev1 * stdDev2 * 0.4)
			return -1.0;
	}
	
	score /= (stdDev1*stdDev2*pixelCount);
	return score;
}

Color Img::getDivergenceAtPixel(int x, int y)
{
	int count = 0;
	Color sum(0.0, 0.0, 0.0);
	
	if(x-1>0)
	{
		sum += pix[(x-1)+y*w];
		count += 1;
	}
	if(x+1<w)
	{
		sum += pix[(x+1)+y*w];
		count += 1;
	}
	if(y-1>0)
	{
		sum += pix[x+(y-1)*w];
		count += 1;
	}
	if(y+1<h)
	{
		sum += pix[x+(y+1)*w];
		count += 1;
	}
	
	sum -= pix[x+y*w]*count;
	
	return sum;
}
