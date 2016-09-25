/*
 *  LineDetector.cpp
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 7/13/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "LineDetector.h"
#include <cmath>
#include <fstream>

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define MAX_JUMP_PIXELS 2
#define START_LINE_LENGTH 3
#define HUGE 999999999.0

void LineDetector::read2DLineSegments(const char *filename, vector<Line> &lineSegments, float thresh)
{
	ifstream fin(filename, ios::in);
	Vec2f start, end;
	int noLines;
	fin >> noLines;
	
	for(int i=0; i<noLines; i++)
	{
		fin >> start[0] >> start[1] >> end[0] >> end[1];
		if((start-end).length() > thresh)
		{
			Line l(start, end);
			lineSegments.push_back(l);
		}
	}
	fin.close();
}

void LineDetector::save2DLineSegments(const char *filename, vector<Line> &lineSegments, float thresh)
{
	ofstream fout(filename, ios::out);
	Vec2f start, end;
	fout << lineSegments.size() << endl;
	for(int i=0; i<lineSegments.size(); i++)
	{
		start = lineSegments[i].getStartPoint2D();
		end = lineSegments[i].getEndPoint2D();
		if((start-end).length() > thresh)
		{
			fout << start[0] << " " << start[1] << " " << end[0] << " " << end[1] << endl;
		}
	}
	fout.close();
}

void LineDetector::find2DLineSegments(Photo *p, float sigma, float tlow, float thigh, float lengthThresh, vector<Line> &lineSegments, const char *filename)
{
	unsigned char *edges;
	vector<vector<int> > xPixels;
	vector<vector<int> > yPixels;
	vector<vector<Line> > linkedLineSegments;
	
	int width = p->getWidth(0);
	int height = p->getHeight(0);
	
	cannyEdgeDetection(p, sigma, tlow, thigh, &edges);
	
	linkEdgePixels(width, height, edges, 50, xPixels, yPixels);
	
	Img temp(*(p->getImg(0)));
	segmentLinkedLists(width, height, 2.0, 20, xPixels, yPixels, linkedLineSegments, &temp);
	temp.write("/Users/ceylan/Desktop/lineFinal.png");
	
	for(int i=0; i<linkedLineSegments.size(); i++)
	{
		for(int j=0; j<linkedLineSegments[i].size(); j++)
		{
			lineSegments.push_back(linkedLineSegments[i][j]);
		}
	}
	
	save2DLineSegments(filename, lineSegments, lengthThresh);
}

void LineDetector::findLinesIntersectingRegion(vector<Line> &edges2D, vector<Line> &lines, Vec2f center, float regionWidth, float regionHeight)
{
	int noLines = edges2D.size();
	
	//write the line equations for the edges of the boundary region
	Vec3f line1, line2, line3, line4;
	Vec2f n;
	n[0] = regionWidth; n[1] = 0.0;
	n.normalize();
	line1[0] = n[0]; line1[1] = n[1];
	line1[2] = 0.0 - n[0]*(center[0]-regionWidth/2.0) - n[1]*(center[1]-regionHeight/2.0);
	
	n[0] = 0.0; n[1] = regionHeight;
	n.normalize();
	line2[0] = n[0]; line2[1] = n[1];
	line2[2] = 0.0 - n[0]*(center[0]-regionWidth/2.0) - n[1]*(center[1]-regionHeight/2.0);
	
	n[0] = -regionWidth; n[1] = 0.0;
	n.normalize();
	line3[0] = n[0]; line3[1] = n[1];
	line3[2] = 0.0 - n[0]*(center[0]+regionWidth/2.0) - n[1]*(center[1]-regionHeight/2.0);
	
	n[0] = 0.0; n[1] = -regionHeight;
	n.normalize();
	line4[0] = n[0]; line4[1] = n[1];
	line4[2] = 0.0 - n[0]*(center[0]+regionWidth/2.0) - n[1]*(center[1]+regionHeight/2.0);
	
	vector<Line*> correspondingLines;
	
	//for each region find the intersecting lines
	Vec2f startPt, endPt;
	
	for(int j=0; j<noLines; j++)
	{
		startPt = edges2D[j].getStartPoint2D();
		endPt = edges2D[j].getEndPoint2D();
		
		bool intersect1 = false;
		bool intersect2 = false;
		bool intersect3 = false;
		bool intersect4 = false;
		Vec2f p;
		float t = (0.0 - line1[2] - dot(Vec2f(line1[0], line1[1]), startPt)) / 
					dot(Vec2f(line1[0], line1[1]), (endPt-startPt));
		if(t>=0.0 && t<=1.0)
		{
			p = startPt + (endPt-startPt)*t;
			if(p[1]>=center[1]-regionHeight/2.0 && p[1]<=center[1]+regionHeight/2.0)
			{
				intersect1 = true;
				//clip
				float distCenter = dot(Vec2f(line1[0], line1[1]), center) + line1[2];
				float distStart = dot(Vec2f(line1[0], line1[1]), startPt) + line1[2];
				float distEnd = dot(Vec2f(line1[0], line1[1]), endPt) + line1[2];
				if(distStart*distCenter < 0.0)
				{
					edges2D[j].setClippedStartPoint(startPt + t*(endPt-startPt));
					edges2D[j].setClippedEndPoint(endPt);
				}
				else
				{
					edges2D[j].setClippedEndPoint(startPt + t*(endPt-startPt));
					edges2D[j].setClippedStartPoint(startPt);
				}
			}
		}
		
		t = (0.0 - line2[2] - dot(Vec2f(line2[0], line2[1]), startPt)) / 
				dot(Vec2f(line2[0], line2[1]), (endPt-startPt));
		if(t>=0.0 && t<=1.0)
		{
			p = startPt + (endPt-startPt)*t;
			if(p[0]>=center[0]-regionWidth/2.0 && p[0]<=center[0]+regionWidth/2.0)
			{
				intersect2 = true;
				//clip
				float distCenter = dot(Vec2f(line2[0], line2[1]), center) + line2[2];
				float distStart = dot(Vec2f(line2[0], line2[1]), startPt) + line2[2];
				float distEnd = dot(Vec2f(line2[0], line2[1]), endPt) + line2[2];
				if(distStart*distCenter < 0.0)
				{
					edges2D[j].setClippedStartPoint(startPt + t*(endPt-startPt));
					edges2D[j].setClippedEndPoint(endPt);
				}
				else
				{
					edges2D[j].setClippedEndPoint(startPt + t*(endPt-startPt));
					edges2D[j].setClippedStartPoint(startPt);
				}
			}
		}
		
		t = (0.0 - line3[2] - dot(Vec2f(line3[0], line3[1]), startPt)) / 
				dot(Vec2f(line3[0], line3[1]), (endPt-startPt));
		if(t>=0.0 && t<=1.0)
		{
			p = startPt + (endPt-startPt)*t;
			if(p[1]>=center[1]-regionHeight/2.0 && p[1]<=center[1]+regionHeight/2.0)
			{
				intersect3 = true;
				//clip
				float distCenter = dot(Vec2f(line3[0], line3[1]), center) + line3[2];
				float distStart = dot(Vec2f(line3[0], line3[1]), startPt) + line3[2];
				float distEnd = dot(Vec2f(line3[0], line3[1]), endPt) + line3[2];
				if(distStart*distCenter < 0.0)
				{
					edges2D[j].setClippedStartPoint(startPt + t*(endPt-startPt));
					edges2D[j].setClippedEndPoint(endPt);
				}
				else
				{
					edges2D[j].setClippedEndPoint(startPt + t*(endPt-startPt));
					edges2D[j].setClippedStartPoint(startPt);
				}
			}
		}
		
		t = (0.0 - line4[2] - dot(Vec2f(line4[0], line4[1]), startPt)) / 
				dot(Vec2f(line4[0], line4[1]), (endPt-startPt));
		if(t>=0.0 && t<=1.0)
		{
			p = startPt + (endPt-startPt)*t;
			if(p[0]>=center[0]-regionWidth/2.0 && p[0]<=center[0]+regionWidth/2.0)
			{
				intersect4 = true;
				//clip
				float distCenter = dot(Vec2f(line4[0], line4[1]), center) + line4[2];
				float distStart = dot(Vec2f(line4[0], line4[1]), startPt) + line4[2];
				float distEnd = dot(Vec2f(line4[0], line4[1]), endPt) + line4[2];
				if(distStart*distCenter < 0.0)
				{
					edges2D[j].setClippedStartPoint(startPt + t*(endPt-startPt));
					edges2D[j].setClippedEndPoint(endPt);
				}
				else
				{
					edges2D[j].setClippedEndPoint(startPt + t*(endPt-startPt));
					edges2D[j].setClippedStartPoint(startPt);
				}
			}
		}
		
		if(intersect1 || intersect2 || intersect3 || intersect4)
		{
			
			lines.push_back(edges2D[j]);
		}
		else 
		{
			//check if line is totally inside
			if(startPt[0]>=center[0]-regionWidth/2.0 && startPt[0]<=center[0]+regionWidth/2.0
			   && startPt[1]>=center[1]-regionHeight/2.0 && startPt[1]<=center[1]+regionHeight/2.0
			   && endPt[0]>=center[0]-regionWidth/2.0 && endPt[0]<=center[0]+regionWidth/2.0
			   && endPt[1]>=center[1]-regionHeight/2.0 && endPt[1]<=center[1]+regionHeight/2.0)
			{
				edges2D[j].setClippedStartPoint(startPt);
				edges2D[j].setClippedEndPoint(endPt);
				lines.push_back(edges2D[j]);
				
			}
		}
	}
}

int LineDetector::findAvgCompatibleEdge(vector<Line> &lines, Vec2f startPt, Vec2f endPt, Vec3f &avgLine, vector<int> &closeEdges, 
										int boundaryThresh, Vec2f translation)
{
	int noLines = lines.size();
	float angleThresh = 5.0*M_PI/180.0;
	
	Vec2f direction = endPt - startPt;
	int length = direction.length();
	direction.normalize();
	
	Vec2f normal((endPt-startPt)[1], (startPt-endPt)[0]);
	normal.normalize();
	Vec3f lineEq;
	lineEq[0] = normal[0]; lineEq[1] = normal[1];
	lineEq[2] = 0.0 - normal[0]*(endPt[0]+translation[0]) - normal[1]*(endPt[1]+translation[1]);

	int count = 0;
	Vec2f avgNormal(0.0, 0.0);
	float avgOffset = 0;
	bool refLineSet = false;
	Vec2f refLineNormal;
	
	Vec2f start, end;
	for(int k=0; k<noLines; k++)
	{
		start = lines[k].getStartPoint2D();
		end = lines[k].getEndPoint2D();
		
		Vec2f lN((end-start)[1], (start-end)[0]);
		lN.normalize();
		
		if(abs(dot(normal, lN)) > cos(angleThresh))
		{
			float d = (abs(dot(lineEq, expand2To3(start))) + abs(dot(lineEq, expand2To3(end)))) / 2.0;
			if(d > boundaryThresh)
				continue;
			
			//check if the projection of the line intersects
			float lStart = dot(direction, start - startPt - translation);
			float lEnd = dot(direction, end - startPt - translation);
			
			if((lStart>length && lEnd>length) || (lStart<0.0 && lEnd<0.0))
				continue;
			
			closeEdges.push_back(k);
			if(!refLineSet)
			{
				refLineNormal = lN;
				refLineSet = true;
			}
			else
			{
				if(dot(refLineNormal, lN) < dot(refLineNormal, lN*-1.0))
					lN *= -1.0;
			}
			
			//find line offset
			float offset = 0.0 - dot(lN, start);
			avgNormal += lN;
			avgOffset += offset;
			count++;
		}
	}
	
	if(count != 0)
	{
		avgNormal /= count;
		avgOffset /= count;
		avgLine[0] = avgNormal[0]; avgLine[1] = avgNormal[1]; avgLine[2] = avgOffset;
	}
	return count;
}


int LineDetector::findBestCompatibleEdge(vector<Line> &lines, Vec2f startPt, Vec2f endPt, float &score, int boundaryThresh, Vec2f translation)
{
	int noLines = lines.size();
	float angleThresh = 5.0*M_PI/180.0;
	
	Vec2f direction = endPt - startPt;
	int length = direction.length();
	direction.normalize();
	
	Vec2f normal((endPt-startPt)[1], (startPt-endPt)[0]);
	normal.normalize();
	Vec3f lineEq;
	lineEq[0] = normal[0]; lineEq[1] = normal[1];
	lineEq[2] = 0.0 - normal[0]*(endPt[0]+translation[0]) - normal[1]*(endPt[1]+translation[1]);
	
	int index = -1;
	float dist = 0.0;
	float maxScore = 0.0;
	vector<Vec2f> overlaps;
	
	Vec2f start, end;
	for(int k=0; k<noLines; k++)
	{
		start = lines[k].getStartPoint2D();
		end = lines[k].getEndPoint2D();
		
		Vec2f lN((end-start)[1], (start-end)[0]);
		lN.normalize();
		if(abs(dot(normal, lN)) > cos(angleThresh))
		{
			float d = (abs(dot(lineEq, expand2To3(start))) + abs(dot(lineEq, expand2To3(end)))) / 2.0;
			if(d > boundaryThresh)
				continue;
			
			//check if the projection of the line intersects
			float lStart = dot(direction, start - startPt);
			float lEnd = dot(direction, end - startPt);
				
			float projLength = abs(lStart-lEnd);
				
			float lClippedStart = lStart;
			float lClippedEnd = lEnd;
				
			if((lClippedStart>length && lClippedEnd>length) || (lClippedStart<0.0 && lClippedEnd<0.0))
				continue;
				
			if(lClippedEnd < lClippedStart)
			{
				float temp = lClippedEnd;
				lClippedEnd = lClippedStart;
				lClippedStart = temp;
			}
			if(lClippedStart < 0.0)
				lClippedStart = 0.0;
			if(lClippedEnd > length)
				lClippedEnd = length;
				
			//check overlap
			if(abs(lClippedStart-lClippedEnd) < projLength/5)
				continue;
				
			if(index != -1)
			{
				if(d<dist)
				{
					dist = d;
					index = k;
				}
			}
			else 
			{
				dist = d;
				index = k;
					
			}
				
			//add overlap
			vector<int> intersections;
			int rem=0;
			for(int o=0; o<overlaps.size(); o++)
			{
				if((lClippedStart > overlaps[o][1] && lClippedEnd > overlaps[o][1]) ||
					(lClippedStart < overlaps[o][0] && lClippedEnd < overlaps[o][0]))
					continue;
				intersections.push_back(o);
			}
			for(int o=0; o<intersections.size(); o++)
			{
				if(overlaps[intersections[o]-rem][0] < lClippedStart)
					lClippedStart = overlaps[intersections[o]-rem][0];
				if(overlaps[intersections[o]-rem][1] > lClippedEnd)
					lClippedEnd = overlaps[intersections[o]-rem][1];
				overlaps.erase(overlaps.begin()+intersections[o]-rem);
				rem++;
			}
			overlaps.push_back(Vec2f(lClippedStart, lClippedEnd));
		}
	}
	
	score = 0;
	for(int i=0; i<overlaps.size(); i++)
	{
		score += overlaps[i][1] - overlaps[i][0];
	}
	
	score /= length;
	return index;
}

void LineDetector::findRepeatingPathsWithLineCompatibility(vector<Line> pathLines, vector<Line> &lines, vector<Vec2f> &matches)
{
	int noLines = pathLines.size();
	float minx, maxx, miny, maxy;
	minx = pathLines[0].getStartPoint2D()[0]; 
	maxx = minx;
	miny = pathLines[0].getEndPoint2D()[1];
	maxy = miny;
	vector<float> pathLineLengths;
	float totalLength = 0.0;
	
	for(int i=0; i<noLines; i++)
	{
		if(pathLines[i].getStartPoint2D()[0] < minx)
			minx = pathLines[i].getStartPoint2D()[0];
		else if(pathLines[i].getStartPoint2D()[0] > maxx)
			maxx = pathLines[i].getStartPoint2D()[0];
		if(pathLines[i].getStartPoint2D()[1] < miny)
			miny = pathLines[i].getStartPoint2D()[1];
		else if(pathLines[i].getStartPoint2D()[1] > maxy)
			maxy = pathLines[i].getStartPoint2D()[1];
		
		if(pathLines[i].getEndPoint2D()[0] < minx)
			minx = pathLines[i].getEndPoint2D()[0];
		else if(pathLines[i].getEndPoint2D()[0] > maxx)
			maxx = pathLines[i].getEndPoint2D()[0];
		if(pathLines[i].getEndPoint2D()[1] < miny)
			miny = pathLines[i].getEndPoint2D()[1];
		else if(pathLines[i].getEndPoint2D()[1] > maxy)
			maxy = pathLines[i].getEndPoint2D()[1];
		
		pathLineLengths.push_back((pathLines[i].getStartPoint2D()-pathLines[i].getEndPoint2D()).length());
		totalLength += (pathLines[i].getStartPoint2D()-pathLines[i].getEndPoint2D()).length();
	}
	
	int border = 5;
	
	int width = (maxx - minx)+border;
	int height = (maxy - miny)+border;
	if(width%2 == 0)
		width += 1;
	if(height%2 == 0)
		height += 1;
	
	int searchWidth = width + 50;
	int searchHeight = height + 50;
	
	Vec2f orgCenter((minx+maxx)/2.0, (miny+maxy)/2.0);
	
	vector<float> scores;
	scores.resize(searchWidth*searchHeight, 0.0);
	
	int boundaryThresh = 10;
	
#pragma omp parallel for
	for(int x=border+width/2; x<searchWidth-border-width/2; x++)
	{
		for(int y=border+height/2; y<searchHeight-border-height/2; y++)
		{
			Vec2f c(x,y);
			
			Vec2f translation = c-orgCenter;
			vector<Line> intersectingLines;
			findLinesIntersectingRegion(lines, intersectingLines, c, width, height );
			
			int count = 0;
			float totalScore = 0.0;
			for(int i=0;i<pathLines.size(); i++)
			{
				float score;
				int index = findBestCompatibleEdge(intersectingLines, pathLines[i].getStartPoint2D()+translation, pathLines[i].getEndPoint2D()+translation, score,  boundaryThresh);
				totalScore += score*pathLineLengths[i];
				if(index != -1)
				{
					count++;
				}
			}
			
			//if(count == noLines)
			{
				scores[y*searchWidth+x] = totalScore/totalLength;
				
			}
		}
	}
	
	/*Img *lineMatching = new Img(rectWidth, rectHeight, &(scores[0]), true);
	lineMatching->write("/Users/ceylan/Desktop/lineMatching.png");
	
	
	int windowSize = 21;
	
	int halfWindowSize = windowSize / 2;
	for( int y = 0; y < rectHeight; y++)
	{
		for( int x = 0; x < rectWidth; x++ )
		{
			if(scores[y*rectWidth + x] == 0.0)
				continue;
			
			bool localMaxima = true;
			
			for(int i=-halfWindowSize; i<=halfWindowSize; i++)
			{
				for(int j=-halfWindowSize; j<=halfWindowSize; j++)
				{
					if(x+i>=0 && x+i<rectWidth && y+j>=0 && y+j<rectHeight)
					{
						if(scores[(y+j)*rectWidth+x+i] > scores[y*rectWidth+x])
						{
							localMaxima = false;
							break;
						}
					}
				}
				if(!localMaxima)
					break;
			}
			if(localMaxima)
			{
				matches.push_back(Vec2f(x,y));
				printf("match at %d %d:%f\n", x, y, scores[y*rectWidth+x]);
			}
		}
	}
	//check for intersections
	vector<Vec2f> newMatches;
	
	for(int i=0; i<matches.size(); i++)
	{
		int upperx = matches[i][0]-(width+border)/2;
		int uppery = matches[i][1]-(height+border)/2;
		int lowerx = matches[i][0]+(width+border)/2;
		int lowery = matches[i][1] + (height+border)/2;
		bool add=true;
		for(int j=0; j<matches.size(); j++)
		{
			if(matches[j][0]>upperx && matches[j][0]<lowerx && matches[j][1]>uppery && matches[j][1]<lowery)
			{
				if(scores[matches[j][1]*rectWidth+matches[j][0]] > scores[matches[i][1]*rectWidth+matches[i][0]])
				{
					add = false;
					break;
				}
			}
		}
		if(add)
		{
			if(scores[matches[i][1]*rectWidth+matches[i][0]] > 0.5)
			{
				newMatches.push_back(matches[i]);
				printf("finalMatch at %f %f:%f\n", matches[i][0], matches[i][1], scores[matches[i][1]*rectWidth+matches[i][0]]);
			}
			
		}
	}
	matches.clear();
	matches.assign(newMatches.begin(), newMatches.end());*/
}

void LineDetector::cannyEdgeDetection(Photo *p, float sigma, float tlow, float thigh, unsigned char **edges)
{
	int width = p->getWidth(0);
	int height = p->getHeight(0);
	
	unsigned char *nms;
	short *derX;
	short *derY;
	short *magnitude;
	
	float *smoothedImgR;
	float *smoothedImgG;
	float *smoothedImgB;
	
	(*edges) = new unsigned char[width*height];
	nms = new unsigned char[width*height];
	
	//perform gaussian smoothing
	(p->getImg(0))->gaussianSmooth(sigma, &smoothedImgR, &smoothedImgG, &smoothedImgB);
	
	//compute gradients
	computeDerivatives(smoothedImgR, smoothedImgG, smoothedImgB, width, height, &derX, &derY);
	
	//compute gradient magnitudes
	computeGradientMagnitudes(derX, derY, width, height, &magnitude);
	
	//compute non-max-suppresion
	computeNonMaxSuppresion(magnitude, derX, derY, width, height, nms);
	
	//apply hysteresis
	applyHysteresis(magnitude, nms, width, height, tlow, thigh, *edges);
	
	delete smoothedImgR;
	delete smoothedImgG;
	delete smoothedImgB;
	delete derX;
	delete derY;
	delete magnitude;
	delete nms;	
}

void LineDetector::computeGradientMagnitudes(short *derX, short *derY, int width, int height, short **magnitude)
{
	int x, y, pos, sq1, sq2;
	
	(*magnitude) = new short [width*height];
	
	for(y=0,pos=0;y<height;y++)
	{
		for(x=0;x<width;x++,pos++)
		{
			sq1 = (int)derX[pos] * (int)derX[pos];
			sq2 = (int)derY[pos] * (int)derY[pos];
			(*magnitude)[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
		}
	}
	
}

void LineDetector::computeDerivatives(float *smoothedImgR, float *smoothedImgG, float *smoothedImgB, int width, int height, short **derX, short **derY)
{
	int x, y, pos;
	
	(*derX) = new short [width*height];
	(*derY) = new short [width*height];
	
	float a, b, c;
	
	//x-derivative
	for(y=0;y<height;y++)
	{
		pos = y * width;
		a = smoothedImgR[pos+1] - smoothedImgR[pos];
		b = smoothedImgG[pos+1] - smoothedImgG[pos];
		c = smoothedImgB[pos+1] - smoothedImgB[pos];
		if(a>b && a>c)
			(*derX)[pos] = a*255;
		else if(b>a && c>a)
			(*derX)[pos] = b*255;
		else
			(*derX)[pos] = c*255;
		//(*derX)[pos] = smoothedImg[pos+1] - smoothedImg[pos]) * 255;
		pos++;
		for(x=1;x<(width-1);x++,pos++)
		{
			a = smoothedImgR[pos+1] - smoothedImgR[pos-1];
			b = smoothedImgG[pos+1] - smoothedImgG[pos-1];
			c = smoothedImgB[pos+1] - smoothedImgB[pos-1];
			if(a>b && a>c)
				(*derX)[pos] = a*255;
			else if(b>a && c>a)
				(*derX)[pos] = b*255;
			else
				(*derX)[pos] = c*255;
			//(*derX)[pos] = (smoothedImg[pos+1] - smoothedImg[pos-1]) * 255;
		}
		a = smoothedImgR[pos] - smoothedImgR[pos-1];
		b = smoothedImgG[pos] - smoothedImgG[pos-1];
		c = smoothedImgB[pos] - smoothedImgB[pos-1];
		if(a>b && a>c)
			(*derX)[pos] = a*255;
		else if(b>a && c>a)
			(*derX)[pos] = b*255;
		else
			(*derX)[pos] = c*255;
		//(*derX)[pos] = (smoothedImg[pos] - smoothedImg[pos-1]) * 255;
	}
	
	//y-derivative
	for(x=0;x<width;x++)
	{
		pos = x;
		a = smoothedImgR[pos+width] - smoothedImgR[pos];
		b = smoothedImgG[pos+width] - smoothedImgG[pos];
		c = smoothedImgB[pos+width] - smoothedImgB[pos];
		if(a>b && a>c)
			(*derY)[pos] = a*255;
		else if(b>a && c>a)
			(*derY)[pos] = b*255;
		else
			(*derY)[pos] = c*255;
		//(*derY)[pos] = (smoothedImg[pos+width] - smoothedImg[pos])*255;
		pos += width;
		for(y=1;y<(height-1);y++,pos+=width)
		{
			a = smoothedImgR[pos+width] - smoothedImgR[pos-width];
			b = smoothedImgG[pos+width] - smoothedImgG[pos-width];
			c = smoothedImgB[pos+width] - smoothedImgB[pos-width];
			if(a>b && a>c)
				(*derY)[pos] = a*255;
			else if(b>a && c>a)
				(*derY)[pos] = b*255;
			else
				(*derY)[pos] = c*255;
			//(*derY)[pos] = (smoothedImg[pos+width] - smoothedImg[pos-width])*255;
		}
		a = smoothedImgR[pos] - smoothedImgR[pos-width];
		b = smoothedImgG[pos] - smoothedImgG[pos-width];
		c = smoothedImgB[pos] - smoothedImgB[pos-width];
		if(a>b && a>c)
			(*derY)[pos] = a*255;
		else if(b>a && c>a)
			(*derY)[pos] = b*255;
		else
			(*derY)[pos] = c*255;
		//(*derY)[pos] = (smoothedImg[pos] - smoothedImg[pos-width])*255;
	}
}

void LineDetector::followEdges(unsigned char *edgemapptr, short *edgemagptr, short lowval, int width)
{
	short *tempmagptr;
	unsigned char *tempmapptr;
	int i;
	int x[8] = {1,1,0,-1,-1,-1,0,1}, y[8] = {0,1,1,1,0,-1,-1,-1};
	
	for(i=0;i<8;i++){
		tempmapptr = edgemapptr - y[i]*width + x[i];
		tempmagptr = edgemagptr - y[i]*width + x[i];
		
		if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
			*tempmapptr = (unsigned char) EDGE;
			followEdges(tempmapptr,tempmagptr, lowval, width);
		}
	}
}

void LineDetector::applyHysteresis(short *mag, unsigned char *nms, int width, int height, float tlow, float thigh, unsigned char *edge)
{
	Img *temp = new Img(width, height);
	
	int x, y, pos, numedges, highcount, lowthreshold, highthreshold;
	int hist[32768];
	short maximum_mag;
	
	for(y=0,pos=0;y<height;y++)
	{
		for(x=0;x<width;x++,pos++)
		{
			if(nms[pos] == POSSIBLE_EDGE) 
				edge[pos] = POSSIBLE_EDGE;
			else edge[pos] = NOEDGE;
		}
	}
	
	//boundary
	for(y=0,pos=0;y<height;y++,pos+=width)
	{
		edge[pos] = NOEDGE;
		edge[pos+width-1] = NOEDGE;
	}
	pos = (height-1) * width;
	for(x=0;x<width;x++,pos++){
		edge[x] = NOEDGE;
		edge[pos] = NOEDGE;
	}
	
	//histogram
	for(int r=0;r<32768;r++) hist[r] = 0;
	for(y=0,pos=0;y<height;y++)
	{
		for(x=0;x<width;x++,pos++)
		{
			if(edge[pos] == POSSIBLE_EDGE) 
				hist[mag[pos]]++;
		}
	}
	
	numedges = 0;
	for(int r=1;r<32768;r++)
	{
		if(hist[r] != 0) 
			maximum_mag = r;
		numedges += hist[r];
	}
	
	highcount = (int)(numedges * thigh + 0.5);
	
	int r = 1;
	numedges = hist[1];
	while((r<(maximum_mag-1)) && (numedges < highcount)){
		r++;
		numedges += hist[r];
	}
	highthreshold = r;
	lowthreshold = (int)(highthreshold * tlow + 0.5);
	
	for(y=0,pos=0;y<height;y++){
		for(x=0;x<width;x++,pos++){
			if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
				edge[pos] = EDGE;
				followEdges((edge+pos), (mag+pos), lowthreshold, width);
			}
		}
	}
	
	for(y=0,pos=0;y<height;y++)
	{
		for(x=0;x<width;x++,pos++) 
		{
			if(edge[pos] != EDGE) 
			{
				edge[pos] = NOEDGE;
				temp->setColor(x,y,Color(0.0,0.0,0.0));
			}
			else
				temp->setColor(x,y,Color(1.0,1.0,1.0));
		}
	}
	//temp->write("/Users/ceylan/Desktop/lines.png");
}

void LineDetector::computeNonMaxSuppresion(short *magnitude, short *derX, short *derY, int width, int height, unsigned char *result)
{
    int rowcount, colcount,count;
    short *magrowptr,*magptr;
    short *gxrowptr,*gxptr;
    short *gyrowptr,*gyptr,z1,z2;
    short m00,gx,gy;
    float mag1,mag2,xperp,yperp;
    unsigned char *resultrowptr, *resultptr;
	
	//Zero the first and last rows of the result image.
	for(count=0,resultrowptr=result,resultptr=result+width*(height-1); 
        count<width; resultptr++,resultrowptr++,count++)
	{
        *resultrowptr = *resultptr = (unsigned char) 0;
    }
	
	//Zero the first and last columns of the result image
    for(count=0,resultptr=result,resultrowptr=result+width-1;
        count<height; count++,resultptr+=width,resultrowptr+=width)
	{
        *resultptr = *resultrowptr = (unsigned char) 0;
    }
	
	//Suppress non-maximum points.
	for(rowcount=1,magrowptr=magnitude+width+1,gxrowptr=derX+width+1,gyrowptr=derY+width+1,resultrowptr=result+width+1;
		rowcount<height-2; 
		rowcount++,magrowptr+=width,gyrowptr+=width,gxrowptr+=width,resultrowptr+=width)
	{   
		for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,resultptr=resultrowptr;
			colcount<width-2; 
			colcount++,magptr++,gxptr++,gyptr++,resultptr++)
		{   
			m00 = *magptr;
			if(m00 == 0)
			{
				*resultptr = (unsigned char) NOEDGE;
			}
			else
			{
				xperp = -(gx = *gxptr)/((float)m00);
				yperp = (gy = *gyptr)/((float)m00);
			}
			
			if(gx >= 0)
			{
				if(gy >= 0)
				{
                    if (gx >= gy)
                    {  
                        /* 111 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - width - 1);
						
                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                        
                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + width + 1);
						
                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {    
                        /* 110 */
                        /* Left point */
                        z1 = *(magptr - width);
                        z2 = *(magptr - width - 1);
						
                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;
						
                        /* Right point */
                        z1 = *(magptr + width);
                        z2 = *(magptr + width + 1);
						
                        mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp; 
                    }
                }
                else
                {
                    if (gx >= -gy)
                    {
                        /* 101 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + width - 1);
						
                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;
						
                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - width + 1);
						
                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {    
                        /* 100 */
                        /* Left point */
                        z1 = *(magptr + width);
                        z2 = *(magptr + width - 1);
						
                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;
						
                        /* Right point */
                        z1 = *(magptr - width);
                        z2 = *(magptr - width + 1);
						
                        mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp; 
                    }
                }
            }
            else
            {
                if ((gy = *gyptr) >= 0)
                {
                    if (-gx >= gy)
                    {          
                        /* 011 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - width + 1);
						
                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;
						
                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + width - 1);
						
                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /* 010 */
                        /* Left point */
                        z1 = *(magptr - width);
                        z2 = *(magptr - width + 1);
						
                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;
						
                        /* Right point */
                        z1 = *(magptr + width);
                        z2 = *(magptr + width - 1);
						
                        mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                    }
                }
                else
                {
                    if (-gx > -gy)
                    {
                        /* 001 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + width + 1);
						
                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;
						
                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - width - 1);
						
                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /* 000 */
                        /* Left point */
                        z1 = *(magptr + width);
                        z2 = *(magptr + width + 1);
						
                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;
						
                        /* Right point */
                        z1 = *(magptr - width);
                        z2 = *(magptr - width - 1);
						
                        mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                    }
                }
            } 
			
            /* Now determine if the current point is a maximum point */
			
            if ((mag1 > 0.0) || (mag2 > 0.0))
            {
                *resultptr = (unsigned char) NOEDGE;
            }
            else
            {    
                if (mag2 == 0.0)
				{
                    *resultptr = (unsigned char) NOEDGE;
				}
                else
				{
                    *resultptr = (unsigned char) POSSIBLE_EDGE;
				}
            }
        } 
    }
}

void LineDetector::cleanEdgeImage(unsigned char *image, int width, int height)
{
	int loop1,loop2;
	unsigned char i1,i2,i3,i4;
	
	for(loop1=1;loop1<height-1;loop1++)
	{
        for(loop2=1;loop2<width-1;loop2++)
		{
            if(image[loop1*width + loop2] == EDGE)
			{
				i1 = image[(loop1-1)*width + loop2];
				i2 = image[loop1*width + loop2-1];
				i3 = image[(loop1+1)*width + loop2];
				i4 = image[loop1*width + loop2+1];
				if((i1 == EDGE) && (i2 == EDGE))
				{
             	    image[loop1*width+loop2] = NOEDGE;
				}
				else if((i2 == EDGE) && (i3 == EDGE))
				{
					image[loop1*width+loop2] = NOEDGE;
				}
				else if((i3 == EDGE) && (i4 == EDGE))
				{
             	    image[loop1*width+loop2] = NOEDGE;
				}
				else if((i4 == EDGE) && (i1 == EDGE))
				{
             	    image[loop1*width+loop2] = NOEDGE;
				}
			}
		}
	}
}

void LineDetector::removeIsolatedPoints(unsigned char *edges, int width, int height)
{
	int loop1,loop2;
	unsigned char i1,i2,i3,i4,i6,i7,i8,i9;
	
	for(loop1=1;loop1<height-1;loop1++)
	{
		for(loop2=1;loop2<width-1;loop2++)
		{
      	    if(edges[loop1*width+loop2] == EDGE)
			{
				i1 = edges[(loop1-1)*width+loop2-1];
				i2 = edges[loop1*width+loop2-1];
				i3 = edges[(loop1+1)*width+loop2-1];
				i4 = edges[(loop1-1)*width+loop2];
				i6 = edges[(loop1+1)*width+loop2];
				i7 = edges[(loop1-1)*width+loop2+1];
				i8 = edges[loop1*width+loop2+1];
				i9 = edges[(loop1+1)*width+loop2+1];
				if((i1+i2+i3+i4+i6+i7+i8+i9) == (8*NOEDGE))
				{
             	    edges[loop1*width+loop2] = NOEDGE;
             	}
            } 
		}
	}
}

void LineDetector::linkOpenEdges(unsigned char *edges, int width, int height, int minLength, vector<vector<int> > &xPixels, vector<vector<int> > &yPixels, Img *im)
{
    int loop1,loop2,loop3;
    unsigned char i1,i2,i3,i4,i6,i7,i8,i9;
    int xp,yp;
    bool endOfLine;
	int weight;
	vector<int> xpix, ypix;
	
    for(loop1=0;loop1<height;loop1++)
	{
		for(loop2=0;loop2<width;loop2++)
		{
      	    if(edges[loop1*width+loop2] == EDGE)
			{
				i1 = 0;
				i2 = 0;
				i3 = 0;
				i4 = 0;
				i6 = 0;
				i7 = 0;
				i8 = 0;
				i9 = 0;
				if(edges[(loop1-1)*width+loop2-1] == EDGE)
            	    i1 = 1;
				if(edges[loop1*width+loop2-1] == EDGE)
					i2 = 1;
				if(edges[(loop1+1)*width+loop2-1] == EDGE)
					i3 = 1;
				if(edges[(loop1-1)*width+loop2] == EDGE)
					i4 = 1;
				if(edges[(loop1+1)*width+loop2] == EDGE)
					i6 = 1;
				if(edges[(loop1-1)*width+loop2+1] == EDGE)
					i7 = 1;
				if(edges[loop1*width+loop2+1] == EDGE)
					i8 = 1;
				if(edges[(loop1+1)*width+loop2+1] == EDGE)
					i9 = 1;
				
				//one connected edge found
				if((i1+i2+i3+i4+i6+i7+i8+i9) == 1)
				{
					weight = 0;
					endOfLine = false;
					xpix.clear();
					ypix.clear();
					
					/* track to end of line */
					xp = loop2;
					yp = loop1;
					
					do
					{
						weight += 1;
						xpix.push_back(xp);
						ypix.push_back(yp);
						edges[yp*width+xp] = NOEDGE;
						
						/* goto next pixel if an edge pixel */
						i1 = edges[(yp-1)*width+xp-1];
						i2 = edges[yp*width+xp-1];
						i3 = edges[(yp+1)*width+xp-1];
						i4 = edges[(yp-1)*width+xp];
						i6 = edges[(yp+1)*width+xp];
						i7 = edges[(yp-1)*width+xp+1];
						i8 = edges[yp*width+xp+1];
						i9 = edges[(yp+1)*width+xp+1];
						if(i1 == EDGE)
						{
							xp--;
							yp--;
						}
						else if(i2 == EDGE)
						{
							xp--;
						}
						else if(i3 == EDGE)
						{
							yp++;
							xp--;
						}
						else if(i4 == EDGE)
						{
							yp--;
						}
						else if(i6 == EDGE)
						{
							yp++;
						}
						else if(i7 == EDGE)
						{
							yp--;
							xp++;
						}
						else if(i8 == EDGE)
						{
							xp++;
						}
						else if(i9 == EDGE)
						{
							xp++;
							yp++;
						}
						else
							endOfLine = true;
					}while(!endOfLine);
					
					if(weight > minLength)
					{
						xPixels.push_back(xpix);
						yPixels.push_back(ypix);
						
						for(int i=0; i<xpix.size(); i++)
							im->setColor(xpix[i], ypix[i], Color(1.0, 0.0, 0.0));
					}
				}
				
			}
		}
	}
}

void LineDetector::linkClosedEdges(unsigned char *edges, int width, int height, int minLength, vector<vector<int> > &xPixels, vector<vector<int> > &yPixels, Img *im)
{
    int loop1,loop2;
    unsigned char i1,i2,i3,i4,i6,i7,i8,i9;
    int xp,yp;
    bool endOfLine;
	int weight;
	vector<int> xpix, ypix;
	
    for(loop1=0;loop1<height;loop1++)
	{
		for(loop2=0;loop2<width;loop2++)
		{
      	    if(edges[loop1*width+loop2] == EDGE)
			{
				i1 = 0;
				i2 = 0;
				i3 = 0;
				i4 = 0;
				i6 = 0;
				i7 = 0;
				i8 = 0;
				i9 = 0;
				if(edges[(loop1-1)*width+loop2-1] == EDGE)
            	    i1 = 1;
				if(edges[loop1*width+loop2-1] == EDGE)
					i2 = 1;
				if(edges[(loop1+1)*width+loop2-1] == EDGE)
					i3 = 1;
				if(edges[(loop1-1)*width+loop2] == EDGE)
					i4 = 1;
				if(edges[(loop1+1)*width+loop2] == EDGE)
					i6 = 1;
				if(edges[(loop1-1)*width+loop2+1] == EDGE)
					i7 = 1;
				if(edges[loop1*width+loop2+1] == EDGE)
					i8 = 1;
				if(edges[(loop1+1)*width+loop2+1] == EDGE)
					i9 = 1;
				
				//one connected edge found
				if((i1+i2+i3+i4+i6+i7+i8+i9) == 2)
				{					
					weight = 0;
					endOfLine = false;
					xpix.clear();
					ypix.clear();
					
					/* track to end of line */
					xp = loop2;
					yp = loop1;
					
					do
					{
						weight += 1;
						xpix.push_back(xp);
						ypix.push_back(yp);
						edges[yp*width+xp] = NOEDGE;
						
						/* goto next pixel if an edge pixel */
						i1 = edges[(yp-1)*width+xp-1];
						i2 = edges[yp*width+xp-1];
						i3 = edges[(yp+1)*width+xp-1];
						i4 = edges[(yp-1)*width+xp];
						i6 = edges[(yp+1)*width+xp];
						i7 = edges[(yp-1)*width+xp+1];
						i8 = edges[yp*width+xp+1];
						i9 = edges[(yp+1)*width+xp+1];
						if(i1 == EDGE)
						{
							xp--;
							yp--;
						}
						else if(i2 == EDGE)
						{
							xp--;
						}
						else if(i3 == EDGE)
						{
							yp++;
							xp--;
						}
						else if(i4 == EDGE)
						{
							yp--;
						}
						else if(i6 == EDGE)
						{
							yp++;
						}
						else if(i7 == EDGE)
						{
							yp--;
							xp++;
						}
						else if(i8 == EDGE)
						{
							xp++;
						}
						else if(i9 == EDGE)
						{
							xp++;
							yp++;
						}
						else
							endOfLine = true;
					}while(!endOfLine);
					
					if(weight > minLength)
					{
						xPixels.push_back(xpix);
						yPixels.push_back(ypix);
						
						for(int i=0; i<xpix.size(); i++)
							im->setColor(xpix[i], ypix[i], Color(1.0, 0.0, 0.0));
					}
				}
				
			}
		}
	}
}

void LineDetector::linkEdgePixels(int width, int height, unsigned char *edges, int minLength, vector<vector<int> > &xPixels, vector<vector<int> > &yPixels)
{
	Img *temp = new Img(width, height);
	
    /* clean up image, make all object 8 connected */
    cleanEdgeImage(edges, width, height);
    
    /* remove isolated points - of no interest to us here */
    removeIsolatedPoints(edges, width, height);
	
    /* link open edges */
    linkOpenEdges(edges, width, height, minLength, xPixels, yPixels, temp);
    
    /* remove isolated points - of no interest to us here */
    removeIsolatedPoints(edges, width, height);
    
    /* link closed edges */
    linkClosedEdges(edges, width, height, minLength, xPixels, yPixels, temp);
    
    /* remove isolated points - of no interest to us here */
    removeIsolatedPoints(edges, width, height);
	
	//temp->write("/Users/ceylan/Desktop/denemeFinal.png");
}

void LineDetector::segmentLinkedLists(int width, int height, float lineDevThresh, int minLineLength, vector<vector<int> > &xPixels, vector<vector<int> > &yPixels, vector<vector<Line> > &lineSegments, Img *img)
{
	vector<Line> linkedLineSegments;
	
	int nolinkedLists = xPixels.size();
	int first, last;
	int index;
	float dev;
	
	for(int i=0; i<nolinkedLists; i++)
	{
		linkedLineSegments.clear();
		int linkedListLength = xPixels[i].size();
		first = 0;
		last = linkedListLength - 1;
		while(first <= last-minLineLength)
		{
			dev = measureLineDeviation(xPixels[i], yPixels[i], first, last, index);
			while(dev > lineDevThresh)
			{
				last = index;
				dev = measureLineDeviation(xPixels[i], yPixels[i], first, last, index);
			}
			if(last-first >= minLineLength)
			{
				Line newLine(Vec2f(xPixels[i][first], yPixels[i][first]), Vec2f(xPixels[i][last], yPixels[i][last]));
				linkedLineSegments.push_back(newLine);
				for(int j=first; j<last; j++)
				{
					img->setColor(xPixels[i][j], yPixels[i][j], Color(1.0, 0.0, 0.0));
				}
			}
			first = last;
			last = linkedListLength - 1;
		}
		lineSegments.push_back(linkedLineSegments);
	}
}

float LineDetector::measureLineDeviation(vector<int> &xPixels, vector<int> &yPixels, int first, int last, int &index)
{
	float epsilon = 0.00000001;
	float dist = sqrt((xPixels[first]-xPixels[last])*(xPixels[first]-xPixels[last]) + (yPixels[first]-yPixels[last])*(yPixels[first]-yPixels[last]));
	float y1y2, x2x1, c, maxDist, temp;
	
	maxDist = 0.0;
	index = first;
	if(dist > epsilon)
	{
		//Eqn of line joining end pts (x1 y1) and (x2 y2) can be parameterised by
		//x*(y1-y2) + y*(x2-x1) + y2*x1 - y1*x2 = 0
		//(See Jain, Rangachar and Schunck, "Machine Vision", McGraw-Hill 1996. pp 194-196)
		y1y2 = yPixels[first] - yPixels[last];
		x2x1 = xPixels[last] - xPixels[first];
		c = yPixels[last]*xPixels[first] - yPixels[first]*xPixels[last];
		for(int i=first; i<=last; i++)
		{
			temp = abs(xPixels[i]*y1y2 + yPixels[i]*x2x1 + c)/dist;
			if(temp > maxDist)
			{
				maxDist = temp;
				index = i;
			}
		}
	}
	else 
	{
		for(int i=first; i<=last; i++)
		{
			temp = sqrt((xPixels[i] - xPixels[first])*(xPixels[i] - xPixels[first]) + (yPixels[i] - yPixels[first])*(yPixels[i] - yPixels[first]));
			if(temp > maxDist)
			{
				maxDist = temp;
				index = i;
			}
		}
	}
	
	return maxDist;
}