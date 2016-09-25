/*
 *  ImageRectification.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 12/6/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <cfloat>
#include "ImageRectification.h"
#include "EdgeDetector.h"
#include "MathUtils.h"

extern "C" 
void canny(unsigned char *image, int rows, int cols, float sigma,
					  float tlow, float thigh, short int* smoothed, int *edge);


ImageRectification::ImageRectification()
{
	secondHorizontalDirection = false;
}

Matrix3f ImageRectification::rectifyImage(Img **data, Vec2f &vanishingPtDist, bool &secondHorizontalVPFound)
{
	secondHorizontalDirection = false;
	computeEdgeMap();
	processEdgeMap();
	computeDominantOrientations();
	secondHorizontalVPFound = secondHorizontalDirection;
	computeVanishingPoints(vanishingPtDist);
	computeRectifiedImage(1.0, data);
	return H;
}

Matrix3f ImageRectification::rectifyImageWrtSecondHorizontalDirection(Img **data, Vec2f &vanishingPtDist)
{
	orientationDom[0] = orientationDom[2];
	computeVanishingPoints(vanishingPtDist);
	computeRectifiedImage(1.0, data);
	return H;
}

void ImageRectification::computeEdgeMap()
{
	EdgeDetector ed;
	
	width = im->width();
	height = im->height();
	int sz = max(width, height);
	float smoothSigma = 1.5f*max(sz, 1024) / 1024.f;
	
	imageData = new unsigned char[width*height];
	edgeData = new int[width*height];
	smoothData = new short int[width*height];
	
	for(int y=0; y<height; y++)
	{
		for(int x=0; x<width; x++)
		{
			Color col = (*im)(x,y);
			unsigned char a = (unsigned char)(col[0] * 255.0);
			unsigned char b = (unsigned char)(col[1] * 255.0);
			unsigned char c = (unsigned char)(col[2] * 255.0);
			
			imageData[y*width+x] = (19595*a + 38470*b + 7471*c + 32768) >> 16;
		}
	}
	
	//ed.cannyEdgeDetection(im, smoothSigma, 0.3, 0.8, &smoothedR, &smoothedG, &smoothedB, &edgeData);
	
	canny(imageData, height, width, smoothSigma, 0.3, 0.8, smoothData, edgeData);
}

void ImageRectification::processEdgeMap()
{
	const int margin = 10;
		
	edgePoints.clear();
	edgeOrientations.clear();
	edgeMagnitudes.clear();
		
	///step1. compute gradient
	for(int i = margin; i < height - margin; ++i)
	{
		for(int j = margin; j < width - margin; ++j)
		{
			if(edgeData[i*width+j] > 1) 
			{
				//float a = smoothedR[i*width+j+1] - smoothedR[i*width+j-1];
				//float b = smoothedG[i*width+j+1] - smoothedG[i*width+j-1];
				//float c = smoothedB[i*width+j+1] - smoothedB[i*width+j-1];
				
				float dx, dy;
				
				//if(a>b && a>c)
				//	dx = a;
				//else if(b>a && c>a)
				//	dx = b;
				//else
				//	dx = c;
				
				//a = smoothedR[(i+1)*width+j] - smoothedR[(i-1)*width+j];
				//b = smoothedG[(i+1)*width+j] - smoothedG[(i-1)*width+j];
				//c = smoothedB[(i+1)*width+j] - smoothedB[(i-1)*width+j];
				//if(a>b && a>c)
				//	dy = a;
				//else if(b>a && c>a)
				//	dy = b;
				//else
				//	dy = c;
				
				dx = float(smoothData[i*width+j+1] - smoothData[i*width+j-1]);
				dy = float(smoothData[(i+1)*width+j] - smoothData[(i-1)*width+j]);
				
				edgePoints.push_back(j + 0.5f);
				edgePoints.push_back(i + 0.5f);
				float oo = (float) atan2(dy, dx) + M_PI * 0.5f;
				edgeOrientations.push_back(oo);
				float mm = (float) sqrt(dx * dx + dy * dy);
				edgeMagnitudes.push_back(mm);
			}
		}
	}
}

void ImageRectification::computeDominantOrientations()
{
	const float o_threshold = 0.314159265358979323846f * 0.5;
	const int o_num = 36, o_gap = 180 / o_num, o_span = 5;
	
    ///////////////////////////////////////////////
	float ohistmax = 0, filtered_hist = 0;
	int i, edge_point_count = edgeOrientations.size();
	
	orientationHist.clear();
	orientationHist.resize(ORTBINS, 0.0f);
    
    //compute histogram
	for(i = 0; i < edge_point_count; ++i)
	{
		int oint = int(edgeOrientations[i]  * 180 / M_PI);
		oint = ((oint + (o_gap/2) + 360) % 180) / o_gap;
		orientationHist[oint] += edgeMagnitudes[i];
	}
	
    ///////////////////////////////////////////////////
	for(i = 0; i < o_num; ++i)
	{
		float ow = 1.0f / (fabsf(sinf(i * M_PI /o_num)) + fabsf(cosf(i * M_PI / o_num)));
		orientationHist[i] = orientationHist[i]  * ow; 
	}
	
    
    float upper_limit = FLT_MAX, ohist_max[3];
    int dist_ov_min = o_num, j_vertical = 0, vertical[3] = {0, 0, 0};
    orientationDom.clear();
	orientationDom.resize(3);
	
	for(int j = 0; j < 3; ++j)
    {
        ohist_max[j] = 0;
        orientationDom[j] = -1;
        for(int i = 0; i < o_num; ++i)
        {
			float value = orientationHist[i];
            if(value >= upper_limit || value <= ohist_max[j])
				continue;
			
            /////////////////////////////////
            int local_maximum = 1;
            for(int u = -o_span; u <= o_span; ++u)
            {
                if(u!= 0 && orientationHist[(u + i + o_num)%o_num] > value)
                {
                    local_maximum = 0;
                    break;
                }
            }
            if(local_maximum == 0) 
				continue;  
			
            ohist_max[j] = value;
            orientationDom[j] = i;
        }
        if(j == 0)
        {
            upper_limit = 1.0f;
            for(int i = 0;i < o_num; ++i) 
				orientationHist[i] /= ohist_max[0];
            ohist_max[0] = 1.0f;
        }
		else
        {
            upper_limit = ohist_max[j];
        }
        int dist_ov = abs(orientationDom[j] - o_num / 2);
        if(dist_ov < dist_ov_min)
        {
            j_vertical = j;
            dist_ov_min = dist_ov;
            vertical[j]  = 1;
        }
    }
	
	
    //look for vertical directon.
    /////////////
	if(j_vertical != 1)
	{
        int temp = orientationDom[1];
        orientationDom[1] = orientationDom[j_vertical];
        orientationDom[j_vertical] = temp;
        float ftemp = ohist_max[1];
        ohist_max[1] = ohist_max[j_vertical];
        ohist_max[j_vertical] = ftemp;
        //////////////////
        vertical[j_vertical] = vertical[1]; 
        vertical[1] = 1;
	}
	
    printf("\r\n" "Dominant Orientations: %d(%.2f), %d(%.2f), %d(%.2f)\r\n", 
		   orientationDom[0] * o_gap, ohist_max[0],
		   orientationDom[1] * o_gap, ohist_max[1],
		   orientationDom[2] * o_gap, ohist_max[2]);
	
    ////////////////////////////////////////////////////////
    if(vertical[2] == 0 && ohist_max[2] > 0.25f * ohist_max[0])
    {
		secondHorizontalDirection = true;
		printf("Possibly second horizontal VPs\r\n");
    }
	else
    {
        orientationDom[2] = orientationDom[0];
    }
}

void ImageRectification::findSecondHorizontalVP()
{
    const float m_pi =  3.14159265358979323846f;
	const float o_threshold = 0.314159265358979323846f * 0.5;
	const int o_num = 36, o_gap = 180 / o_num, o_span = 2;
	
    ///////////////////////////////////////////////
	float ohistmax = 0, filtered_hist = 0;
	int i, edge_point_count = edgeOrientations.size();
    float ohist[o_num];
	std::fill(ohist, ohist + o_num, 0.0f);
    
    //////////////////////////////////////////
    //compute histogram
	for(i = 0; i < edge_point_count; ++i)
	{
		
        float dvx = edgePoints[2 * i] * vPoints[0].c - vPoints[0].a;
        float dvy = edgePoints[2 * i + 1] * vPoints[0].c - vPoints[0].b;
        float dvq = sqrt(dvx * dvx + dvy * dvy); 
        float dx = cos(edgeOrientations[i]), dy = sin(edgeOrientations[i]);
        if(fabs(dvy * dy + dvx * dx) > 0.995 * dvq)
        {
            filtered_hist ++;//= _edge_magtitudes[i];
        }else
        {
		    int oint = int( edgeOrientations[i]  * 180 / m_pi);
		    oint = ((oint + (o_gap/2) + 360) % 180) / o_gap;
		    ohist[oint] ++;//= _edge_magtitudes[i];
        }
	}
    ///
    float ohist_max = 0;
    int dist_ov_min = o_num;
	
    for(i = 0; i < o_num; ++i)
    {
        if(i == orientationDom[0]) continue;
        int dist_ov = abs(i - o_num / 2);
        if(dist_ov < o_num / 5) continue;
		
        float value = ohist[i];
        if(value <= ohist_max)continue;
        if(value < 0.1 * filtered_hist) continue;
		
        /////////////////////////////////
        int local_maximum = 1;
        float hsum = 0; 
        for(int u = -o_span; u <= o_span; ++u)
        {
            hsum += ohist[(u + i + o_num)%o_num];
            if(u!= 0 && ohist[(u + i + o_num)%o_num] > value)
            {
                local_maximum = 0;
                break;
            }
        }
        if(local_maximum == 0) continue;  
        if(hsum < 0.4 * filtered_hist) continue;
		
        ohist_max = value;
        orientationDom[2] = i;
    }
    printf("Found second Horizontal VP : %d [%f]\r\n",
		   orientationDom[2] * o_gap, ohist_max / filtered_hist);
}

void ImageRectification::computeVanishingPoints(Vec2f &vanishingPtDist)
{
	
	int i, j;
	const float o_threshold = 0.314159265358979323846f * 0.5;
	int edge_point_count = edgeOrientations.size();
	const int o_num = 36, o_gap = 180 / o_num; //from the resolution of initial orientation histogram
	
    float* ohist = &(orientationHist[0]);
	const int blur_window = 1;
	const int blur_weight[blur_window * 2 + 1] = {1, 3, 1};
	
	const int* bweight = blur_weight + blur_window;
	float xc = im->width() * 0.5f;
	float yc = im->height() * 0.5f;
	int filtered_line_range[3] = {0, 0, 0};
	int vanishing_ok[2] = {0, 0};
    vPoints.clear();
	vPoints.resize(2);
    orientationLineCount.clear();
	orientationLineCount.resize(2);
	orientationLineCount[0] = 0;
    orientationLineCount[1] = 0;
    vPoints[0] = Line2D(1.0f, 0.0f, 0.0f);
    vPoints[1] = Line2D(0.0f, 1.0f, 0.0f);
	
    ////////////////////////////////////
	for(i = 0; i < 2; ++i)
	{
		int o_test = orientationDom[i];
		float fo_test = o_test * M_PI / o_num;
		int int_omin, int_omax; 
		for(int_omax = 1; int_omax <= 5; int_omax++)
		{
			if(ohist[(o_test + int_omax) % o_num] < ohist[o_test] * 0.3) 
				break;
		}
		for(int_omin = -1; int_omin >=-5; int_omin--)
		{
			if(ohist[(o_test + int_omin + o_num) % o_num] < ohist[o_test] * 0.3) 
				break;
		}
		float diff_max = (int_omax - 0.5f) * M_PI / o_num;
		float diff_min = (int_omin + 0.5f) * M_PI / o_num;
		
		int bmin, bmax, b_num; 
		///hough transform for the intersect?
		if(i == 0)
		{
			//get all possile range
			//line equation Y - CY = A * (X - CX) + B
			//get histogram of B, and select all local maximums
			bmin = -height;
			bmax = height;
		}
		else
		{
			//horizonatl direction.
			//line equation X - CX = A * (Y  - CY) + B
			bmin = -width;
			bmax =  width ;
		}
		
		b_num = bmax - bmin + 1;
		int sub_o_num = (int_omax - int_omin + 1) * 3;
		vector<double> ldata, ldatax;
		vector<vector<int> > vote2;
		vote2.resize(b_num);
		for(int k=0; k<b_num; k++)
		{
			vote2[k].resize(sub_o_num, 0);
		}
		
		vector<int> bvote(b_num, 0);
		
		for(j = 0; j < edge_point_count; ++j)
		{
			float doj = fmod(edgeOrientations[j] - fo_test + 1.5f * M_PI, M_PI) - 0.5f * M_PI;
			if(doj < diff_min || doj > diff_max) 
				continue;
			float xj = edgePoints[j * 2] - xc;
			float yj = edgePoints[j * 2 + 1] - yc;  
			float b ; 
			if(i == 0)  
				b = yj - tan(edgeOrientations[j]) * xj;
			else        
				b = xj - yj / tan(edgeOrientations[j]);
			
			//////////////////////////////////////////////////
			int int_b = (int) floor(b );
			int sub_oint = int(floor(doj * 3 * o_num / M_PI + 0.5)) - int_omin * 3;
			
			if(int_b >= bmin + 10 && int_b <= bmax - 10)
			{
				for(int k = -blur_window; k <= blur_window ; ++k)
				{
					bvote[int_b - bmin + k] += bweight[k];
					vote2[int_b - bmin + k][sub_oint] += bweight[k];
				}
			}
		}
		
		int bvote_max = 0;
		for(j = 0; j < bvote.size(); ++j)
		{
			bvote_max = max(bvote_max, bvote[j]);
		}
		for(j = 10; j < bvote.size() -10; ++j)
		{
			if( bvote[j] > bvote_max /  8 && 	
			   bvote[j] > bvote[j -1] && bvote[j] > bvote[j + 1]&&  
			   bvote[j] > bvote[j- 2] && bvote[j] > bvote[j + 2]&&  
			   bvote[j] > bvote[j- 3] && bvote[j] > bvote[j + 3]&&  
			   bvote[j] > bvote[j- 4] && bvote[j] > bvote[j + 4]
			   )
			{
				int vote_max = 0, idx_vote_max = 1;
				for(int k = 1; k < sub_o_num -1; ++k)
				{
					if(vote2[j][k] > vote_max)
					{
						vote_max = vote2[j][k];
						idx_vote_max = k;
					}
				}
				vote_max =  vote2[j][idx_vote_max] + (vote2[j][idx_vote_max -1] + vote2[j][idx_vote_max +1]);
				
				//has a dominant orientation
                if(vote_max > bvote[j] / 3)
				{
					float offset = 0.5f * (vote2[j][idx_vote_max + 1] - vote2[j][idx_vote_max -1])/
					(vote_max * 2 -vote2[j][idx_vote_max + 1] - vote2[j][idx_vote_max -1] );
					float ok = fo_test + (idx_vote_max  + offset + 3 * int_omin) * M_PI / (o_num * 3);
					double cok = cos(ok), sok = sin(ok);
					double line[3] = {-sok, cok, 0};
					if(i == 0)  line[2] = - (( j + bmin + 0.5f) * line[1]);
					else        line[2] = - ((j + bmin + 0.5f) * line[0] );
					ldata.insert(ldata.end(), line, line + 3);
				}
			}
			
		}
		
		double vpoint[3];
		int line_count = ldata.size()/3, new_line_count;
		vector<int> lmask(line_count, 0);
        new_line_count = computeVanishingPoint(line_count, (double (*)[3]) (&ldata[0]), vpoint, &lmask[0], 0.02);
		if(new_line_count <= 2) 
			continue;
		
		for(j = 0; j < line_count; ++j)
		{
			if(lmask[j])
			{
				double *pd = &(ldata[j * 3]);
				edgeLines.push_back(Line2D((float)pd[0], float(pd[1]), float(pd[2] - pd[0] * xc - pd[1] * yc)));
			}
		}
		
		vanishing_ok[i] = 1;
	    filtered_line_range[i + 1] = edgeLines.size();
		
	    vpoint[0] += vpoint[2] * xc;
	    vpoint[1] += vpoint[2] * yc;
		
	    if(vpoint[2] < 0) {vpoint[0] *=-1.0; vpoint[1] *= -1.0; vpoint[2] *= -1.0;}
        printf("VP (%.3f, %.3f, %f) (%d out of %d) [%d:%d]\r\n",
			   vpoint[0], vpoint[1], vpoint[2], new_line_count, line_count, int_omin, int_omax);
		
		
        //////////////////////////////////////////////////////////
	    vPoints[i] = Line2D(vpoint);
        orientationLineCount[i] = new_line_count;
		
		
		//get more lines and refine the vanishing point? no...
		
	}//for the two orientations
	
	vanishingPtDist[0] = vPoints[0].a / vPoints[0].c;//sqrt(vpoint[0]*vpoint[0] + vpoint[1]*vpoint[1]);
	vanishingPtDist[1] = vPoints[0].b / vPoints[0].c;
	printf("VP (%.1f, %.1f)(%.1f, %.1f)\r\n", vPoints[0].a / vPoints[0].c, 
		   vPoints[0].b / vPoints[0].c,vPoints[1].a / vPoints[1].c, vPoints[1].b / vPoints[1].c);
    setHomographyFromVP();
}

int  ImageRectification::computeLineIntersection(double line1[3], double line2[3], double pt[3])
{
	pt[0] = line1[1] * line2[2] - line1[2] * line2[1];
	pt[1] = line1[2] * line2[0] - line1[0] * line2[2];
	pt[2] = line1[0] * line2[1] - line1[1] * line2[0];
	double sq = 1.0 / sqrt(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]);
	pt[0] *= sq;
	pt[1] *= sq;
	pt[2] *= sq;
	return 1;
}

int ImageRectification::linesDuplicate(double a[3], double b[3])
{
	return fabs(a[0] * b[0] + a[1] * b[1] + a[2] * b[2])/
			sqrt((a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2])) > 0.999999999;  
}

int ImageRectification::computeVanishingPoint(int num, double line[][3], double vp[], int mask[], double th)
{
	
	//for every two lines, compute their intersection...
	//count the number of lines crossing this point...
	//const double thp = 5; // five pixels if not at infinity?
	
	int count_max = 0;
	//each line can be seen as a 3d point, and their intersection can be seen as a line.....
	
    if(num < 100)
    {
		for(int i = 0; i < num; ++i)
        {
		    for(int j = i + 1; j < num; ++j)
		    {
			    //skip if if two lines are too close..
			    if(linesDuplicate(line[i], line[j])) 
					continue;
				
			    //find the  crossing of  line[i] and line[j]
				
			    double isp[3]; int count = 0;
				if(!computeLineIntersection(line[i], line[j], isp)) 
					continue;
			    for(int k = 0; k < num; ++k)
			    {
				    double dist = fabs(line[k][0] * isp[0] + line[k][1] * isp[1] + line[k][2] * isp[2]);
				    if(dist < th)
				    {
				    	//fcount += exp(-dist / th / 3);
					    count ++;
				    }
			    }
			    if(count > count_max)
			    {
			    	//fcount_max = fcount;
				    vp[0] = isp[0];
				    vp[1] = isp[1]; 
				    vp[2] = isp[2];
					count_max = count;
			    }
		    }
	    }
    }
	else
    {
		double fcount_max = 0.0;
        for(int idx = 0; idx < 1000; ++idx)
        {
	        int i = rand() % num, j = rand() % num;
            if(i == j || linesDuplicate(line[i], line[j])) 
				continue;
			double isp[3], fcount = 0;	 int count = 0;
			if(!computeLineIntersection(line[i], line[j], isp)) 
				continue;
		    for(int k = 0; k < num; ++k)
		    {
			    double dist = fabs(line[k][0] * isp[0] + line[k][1] * isp[1] + line[k][2] * isp[2]);
			    if(dist < th)
			    {
			    	fcount += exp(-dist / th / 3);
				    count ++;
			    }
		    }
		    if(fcount > fcount_max)
		    {
		    	fcount_max = fcount;
			    vp[0] = isp[0];
			    vp[1] = isp[1]; 
			    vp[2] = isp[2];
			    count_max = count;
		    }
            if(idx> 1000 && num * 9 <= 10 * count_max) 
				break;
	    }
    }
	
	//refine_vanishing_point
	if(count_max > 2)
	{
		int idx = 0;
		vector<double> buf(count_max * 3);
		for(int k = 0; k < num && idx < count_max; ++k)
		{
			if(fabs(line[k][0] * vp[0] + line[k][1] * vp[1] + line[k][2] * vp[2]) < th)
			{
				buf[idx] = line[k][0];
				buf[idx + count_max] = line[k][1];
				buf[idx + count_max * 2] = line[k][2];
				idx ++;
				if(mask) mask[k] = 1;
			}
		}
		MathUtils::lastNullVectorT(&buf[0], vp, count_max, 3);
		if(vp[2] < 0) {vp[0] *= -1.0; vp[1] *= -1.0; vp[2] *= - 1.0;}
		return count_max;
	}
	else if(count_max == 2)
	{
        if(mask)
        {
		    for(int k = 0; k < num; ++k)
		    {
			    if(fabs(line[k][0] * vp[0] + line[k][1] * vp[1] + line[k][2] * vp[2]) < th)
			    {
				    mask[k] = 1;
			    }
		    }
        }
		if(vp[2] < 0) {vp[0] *= -1.0; vp[1] *= -1.0; vp[2] *= - 1.0;}
		return 2;
	}
	else
	{
		return 0;
	}
}

void ImageRectification::setHomographyFromVP()
{
    if(vPoints.size() < 2) 
		return;
	
	float xc = width * 0.5f;
	float yc = height * 0.5f;
	
	//set the matrix to identity first.
    H.setIdentity();
	
	////////////////////////////////////////////
    if(fabs(vPoints[0].a) > 1.0f)
    {
	    H[0][1] = vPoints[0].b / vPoints[0].a;	
		H[0][2] = vPoints[0].c / vPoints[0].a;
    }
	else
    {
        double s = vPoints[0].a > 0 ? 1.0 : - 1.0;
	    H[0][0] = vPoints[0].a * s; H[0][1] = vPoints[0].b * s;	H[0][2] = vPoints[0].c * s;
    }
	
    if(fabs(vPoints[1].b) > 1.0f)
    {
        H[1][0] = vPoints[1].a / vPoints[1].b;   
		H[1][2] = vPoints[1].c / vPoints[1].b;
    }else
    {
        double s = vPoints[1].b > 0 ? 1.0 : - 1.0;
        H[1][0] = vPoints[1].a * s; H[1][1] = vPoints[1].b * s;  	H[1][2] = vPoints[1].c * s;
    }
	
	//H[0][0] *= 10.0;
	//H[0][1] *= 10.0;
	//H[0][2] *= 10.0;
	
	//H[1][0] *= 10.0;
	//H[1][1] *= 10.0;
	//H[1][2] *= 10.0;
	
	Vec3f u = H[0];
	Vec3f v = H[1];
	
	double fSq = (u[0]*v[0]-(u[0]*v[2]+u[2]*v[0])*xc+u[2]*v[2]*xc*xc)+(u[1]*v[1]-(u[1]*v[2]+u[2]*v[1])*yc+u[2]*v[2]*yc*yc);
	
	//vp1[3] = {(H[0][0] - xc * H[0][2]) ,	(H[0][1] - yc * H[0][2]) , 	H[0][2]};
	//double vp2[3] = {(H[1][0] - xc * H[1][2]) ,	(H[1][1] - yc * H[1][2]) ,	H[1][2]};
	//double ftq = -(vp1[0] * vp2[0] + vp1[1] * vp2[1]) / vp1[2] / vp2[2];
	
	//fSq /= u[2];
	//fSq /= v[2];
	//double f = sqrt(fSq*-1.0);
	
	//printf("f:%f\n", f);
	//Matrix3f A;
	//A.setIdentity();
	//A[0][0] = f; A[1][1] = f;
	//A[2][0] = xc; A[2][1] = yc;
	
	//Matrix3f M = (A*A.transpose()).getInverseMatrix();
	
	//float k = dot(v, M*v);
	
	//float l = dot(u, M*u);
	
	//float ratio = sqrt(l/k);
	
	//H[0][0] = u[0]*ratio;
	//H[0][1] = u[1]*ratio;
	//H[0][2] = u[2]*ratio;
	
	//H[1][0] = v[0];
	//H[1][1] = v[1];
	//H[1][2] = v[2];
	
	//printf("length u:%f, v:%f ratio:%f\n", H[0].length(), H[1].length(), ratio);
	
	for(int pass = 0; pass < 2; ++pass)
	{
		double vp1[3] = {(H[0][0] - xc * H[0][2]) ,	(H[0][1] - yc * H[0][2]) , 	H[0][2]};
		double vp2[3] = {(H[1][0] - xc * H[1][2]) ,	(H[1][1] - yc * H[1][2]) ,	H[1][2]};
		double ftq = -(vp1[0] * vp2[0] + vp1[1] * vp2[1]) / vp1[2] / vp2[2];
		double F = 2945.495, Fsq = ftq > F * F ? ftq : F * F;
		double n1 = sqrt(vp1[0] * vp1[0]  + vp1[1] * vp1[1]);
		double n2 = sqrt(vp2[0] * vp2[0]  + vp2[1] * vp2[1]);
		double nv1 = sqrt(vp1[0] * vp1[0] / Fsq + vp1[1] * vp1[1] / Fsq + vp1[2] * vp1[2]);
		double nv2 = sqrt(vp2[0] * vp2[0] / Fsq + vp2[1] * vp2[1] / Fsq + vp2[2] * vp2[2]);
		double nvr = nv1 / nv2;
		double temp1 = 1.0 / sqrt(n1 * n2 * nvr);
		double r1 = temp1, r2 = temp1 * nvr;
		//r1 = 1.0; r2 = 1.0;
		H[0][0] *= r1;	H[0][1] *= r1;	H[0][2] *= r1;
		H[1][0] *= r2;	H[1][1] *= r2;	H[1][2] *= r2;
		
		
		H[2][2]  = 1.0  - H[0][2] * xc - H[1][2] * yc;
		
		double  kk  = H[2][2] + H[0][2] * xc + H[1][2] * yc; 
		H[2][0] = kk * xc - (H[0][0] * xc + H[1][0] * yc);
		H[2][1] = kk * yc - (H[0][1] * xc + H[1][1] * yc);
		
		Hinv = H.getInverseMatrix();
		
		if(Hinv[0][0] < 0 || Hinv[1][1] < 0)
		{
			if(Hinv[0][0] < 0)
			{
				H[0][0] = - H[0][0];
				H[0][1] = - H[0][1];	
				H[0][2] = - H[0][2];
			}
			if(Hinv[1][1] < 0)
			{	
				H[1][0] = - H[1][0]; 
				H[1][1] = - H[1][1];	
				H[1][2] = - H[1][2];
			}
		}
		else
		{
			break;
		}
		
	}
	
	printf("H:\n");
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			printf("%f ", H[j][i]);
		}
		printf("\n");
	}
}

void ImageRectification::computeRectifiedImage(float magnification, Img **data)
{
	if(vPoints.size() < 2)
	{
		printf("No vanishing point found for rectification.\n");
	}
	
	int offsetX[2], offsetY[2];
	
	getRectifiedBoundingBox(offsetX, offsetY, &(vPoints[0]));
	//offsetX[0] = 0; offsetX[1] = 0;
	//offsetY[0] = 0; offsetY[1] = 0;
	
	int newWidth = width;
	int newHeight = height;
	
	newWidth += offsetX[1] - offsetX[0];
	newHeight += offsetY[1] - offsetY[0];
	
	if(newWidth < 32 || newHeight < 32)
	{
		printf("Rectified image size too small.\n");
		(*data) = new Img(width, height, Color(0,0,0));
		return ;
	}
	
	adjustRectificationMatrix(-offsetX[0], -offsetY[0]);
	
	if(magnification >= 0.5f)
	{
		printf("magnification factor:%f\n", magnification);
		H[2][0] *= magnification;
		H[2][1] *= magnification;
		H[2][2] *= magnification;
		
		Hinv[0][2] /= magnification;
		Hinv[1][2] /= magnification;
		Hinv[2][2] /= magnification;
		
		newWidth = (int) (newWidth*magnification);
		newHeight = (int) (newHeight*magnification);
	}
	
	(*data) = new Img(newWidth,newHeight);
	generateRectifiedImage(newWidth, newHeight, *data);
}

void ImageRectification::getRectifiedBoundingBox(int offsetX[2], int offsetY[2], Line2D* vps)
{
	double corner[5][2] = {{0, 0}, { width, 0}, {0, height}, {width, height}, {0.5 * width, 0.5 * height}};
	double tcorner[5][3], zz;
	offsetX[0] = offsetX[1] = offsetY[0] = offsetY[1] = 0;
	for(int t = 0; t < 5; ++t)
	{
		zz = corner[t][0] * Hinv[0][2] + corner[t][1] * Hinv[1][2] + Hinv[2][2];
		tcorner[t][2] = zz;
		tcorner[t][0] = (corner[t][0] * Hinv[0][0] + corner[t][1] * Hinv[1][0] + Hinv[2][0]) / zz;
		tcorner[t][1] = (corner[t][0] * Hinv[0][1] + corner[t][1] * Hinv[1][1] + Hinv[2][1]) / zz;
	}
	
	
	float vp_fix_ratio = 0.5f;
	Line2D& vx = vps[0];
	
	if(vx.a / vx.c > 0 && vx.a / vx.c < width)
	{
		float nx = width * 0.5f * (1.0f - vp_fix_ratio)  + vx.a / vx.c * vp_fix_ratio;
		if(tcorner[0][2] < 0 || tcorner[2][2] < 0  )
		{
			corner[0][0] = corner[2][0] = nx;
			for(int t = 0; t < 3; t +=2)
			{
				zz = corner[t][0] * Hinv[0][2] + corner[t][1] * Hinv[1][2] + Hinv[2][2];
				tcorner[t][2] = zz;
				tcorner[t][0] = (corner[t][0] * Hinv[0][0] + corner[t][1] * Hinv[1][0] + Hinv[2][0]) / zz;
				tcorner[t][1] = (corner[t][0] * Hinv[0][1] + corner[t][1] * Hinv[1][1] + Hinv[2][1]) / zz;
			}
		}else if(tcorner[1][2] < 0 || tcorner[3][2] < 0 )
		{
			corner[1][0] = corner[3][0] = nx;
			for(int t = 1; t < 4; t += 2)
			{
				zz = corner[t][0] * Hinv[0][2] + corner[t][1] * Hinv[1][2] + Hinv[2][2];
				tcorner[t][2] = zz;
				tcorner[t][0] = (corner[t][0] * Hinv[0][0] + corner[t][1] * Hinv[1][0] + Hinv[2][0]) / zz;
				tcorner[t][1] = (corner[t][0] * Hinv[0][1] + corner[t][1] * Hinv[1][1] + Hinv[2][1]) / zz;
			}
		}
	}
	
	Line2D& vy = vps[1];
	if(vy.b / vy.c > 0 && vy.b / vy.c < height)
	{
		float ny = height * 0.5f * (1.0f - vp_fix_ratio)   + vy.b / vy.c * vp_fix_ratio;
		if(tcorner[0][2] < 0 || tcorner[1][2] < 0  )
		{
			corner[0][1] = corner[1][1] = ny;
			for(int t = 0; t < 2; t ++)
			{
				zz = corner[t][0] * Hinv[0][2] + corner[t][1] * Hinv[1][2] + Hinv[2][2];
				tcorner[t][2] = zz;
				tcorner[t][0] = (corner[t][0] * Hinv[0][0] + corner[t][1] * Hinv[1][0] + Hinv[2][0]) / zz;
				tcorner[t][1] = (corner[t][0] * Hinv[0][1] + corner[t][1] * Hinv[1][1] + Hinv[2][1]) / zz;
			}
		}else if(tcorner[2][2] < 0 || tcorner[3][2] < 0 )
		{
			corner[2][1] = corner[3][1] = ny;
			for(int t = 2; t < 4; t ++)
			{
				zz = corner[t][0] * Hinv[0][2] + corner[t][1] * Hinv[1][2] + Hinv[2][2];
				tcorner[t][2] = zz;
				tcorner[t][0] = (corner[t][0] * Hinv[0][0] + corner[t][1] * Hinv[1][0] + Hinv[2][0]) / zz;
				tcorner[t][1] = (corner[t][0] * Hinv[0][1] + corner[t][1] * Hinv[1][1] + Hinv[2][1]) / zz;
			}
		}
	}
	
	
	if(tcorner[0][0] < 0 && tcorner[2][0] < 0)
	{
		offsetX[0] = (int)ceil(max(tcorner[0][0], tcorner[2][0]));
	}else if(tcorner[0][0] > 0 && tcorner[2][0] > 0)
	{
		offsetX[0] = (int) floor(min(tcorner[0][0], tcorner[2][0]));
	}
	
	if(tcorner[0][1] < 0 && tcorner[1][1] < 0)
	{
		offsetY[0] = (int)ceil(max(tcorner[0][1], tcorner[1][1]));
	}else if(tcorner[0][1] > 0 && tcorner[1][1] > 0)
	{
		offsetY[0] = (int) floor(min(tcorner[0][1], tcorner[1][1]));
	}
	
	if(tcorner[1][0] < width && tcorner[3][0] < width)
	{
		offsetX[1] = (int)ceil(max(tcorner[1][0], tcorner[3][0])) - width;
	}else if(tcorner[1][0] > width && tcorner[3][0] > width)
	{
		offsetX[1] = (int) floor(min(tcorner[1][0], tcorner[3][0])) - width;
	}
	
	if(tcorner[2][1] < height && tcorner[3][1] < height)
	{
		offsetY[1] = (int)ceil(max(tcorner[2][1], tcorner[3][1])) - height;
	}else if(tcorner[2][1] > height && tcorner[3][1] > height)
	{
		offsetY[1] = (int) floor(min(tcorner[2][1], tcorner[3][1])) - height;
	}
}

void ImageRectification::adjustRectificationMatrix(double dx, double dy)
{
	if(dx == 0.0 && dy == 0.0)
		return;
	
	H[2][0] -= (H[0][0]*dx + H[1][0]*dy);
	H[2][1] -= (H[0][1]*dx + H[1][1]*dy);
	H[2][2] -= (H[0][2]*dx + H[1][2]*dy);
	
	printf("H:\n");
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			printf("%f ", H[j][i]);
		}
		printf("\n");
	}
	
	Hinv = H.getInverseMatrix();
}

void ImageRectification::transformPointToOriginal(float x, float y, float &tx, float &ty)
{
	float tz;
	tz = H[0][2]*x + H[1][2]*y + H[2][2];
	tx = (H[0][0]*x + H[1][0]*y + H[2][0]) / tz;
	ty = (H[0][1]*x + H[1][1]*y + H[2][1]) / tz;
}

void ImageRectification::generateRectifiedImage(int newWidth, int newHeight, Img *data)
{
	#pragma omp parallel for
	for(int y=0; y<newHeight; y++)
	{
		for(int x=0; x<newWidth; x++)
		{
			float u, v;
			transformPointToOriginal(x+0.5f, y+0.5f, u, v);
			if(u>=0 && u<width && v>=0 && v<height)
			{
				data->setColor(x,y,(*im)(u,v));
			}
			else
			{
				data->setColor(x, y, Color(0.0, 0.0, 0.0));
			}
		}
	}
	data->write("/Users/ceylan/Desktop/rect.png");
}