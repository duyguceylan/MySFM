/*
 *  RepetitionFinder.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 12/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include <opencv2/opencv.hpp>
#include "RepetitionFinder.h"
#include "MathUtils.h"

bool lessWrtYCoord(const Vec2f& lhs, const Vec2f& rhs)
{
	return lhs[1] < rhs[1];
}

bool lessWrtCellYCoord(const repetitionCell& lhs, const repetitionCell& rhs)
{
	return lhs.center[1] < rhs.center[1];
}

bool lessWrtXCoord(const Vec2f& lhs, const Vec2f& rhs)
{
	return lhs[0] < rhs[0];
}

bool lessWrtCellXCoord(const repetitionCell& lhs, const repetitionCell& rhs)
{
	return lhs.center[0] < rhs.center[0];
}

Img RepetitionFinder::getAvgTemplate(vector<Img *> &templateImg, vector<repetitionGrid> &matches, vector<Vec2i> &templateSizes, Vec2i orgTempSize)
{
	int noImages = templateImg.size();
	int noMatches = 0;
	int width, height;
	
	//find min width/height
	for(int i=0; i<noImages; i++)
	{
		if(i==0 && templateSizes[i][0] != 0)
		{
			width = templateSizes[i][0];
			height = templateSizes[i][1];
		}
		else if(templateSizes[i][0] < width && templateSizes[i][0] != 0)
		{
			width = templateSizes[i][0];
			height = templateSizes[i][1];
		}
		noMatches += matches[i].gridCells.size();
	}
	
	Img avgImg(width,height);
	
	for(int i = 0; i < noImages; i++)
	{
		float scale = (float)(templateSizes[i][0]) / float(width);
		for(int m=0; m<matches[i].gridCells.size(); m++)
		{
			Img tmp(templateSizes[i][0], templateSizes[i][1]);
			for(int x=0; x<templateSizes[i][0]; x++)
			{
				for(int y=0; y<templateSizes[i][1]; y++)
				{
					if(matches[i].gridCells[m].center[0]-templateSizes[i][0]/2+x >= 0 && matches[i].gridCells[m].center[0]-templateSizes[i][0]/2+x < templateImg[i]->width() &&
					   matches[i].gridCells[m].center[1]-templateSizes[i][1]/2+y >= 0 && matches[i].gridCells[m].center[1]-templateSizes[i][1]/2+y < templateImg[i]->height())
					{
						tmp.setColor(x, y, (*templateImg[i])(matches[i].gridCells[m].center[0]-templateSizes[i][0]/2+x, 
															 matches[i].gridCells[m].center[1]-templateSizes[i][1]/2+y));
					}
					else
					{
						tmp.setColor(x, y, Color(0.0,0.0,0.0));
					}
				}
			}
			
			Img scaledTmp(width, height);
			for(int x=0; x<width; x++)
			{
				for(int y=0; y<height; y++)
				{
					if(x*scale>=0 && x*scale < templateSizes[i][0] && y*scale>=0 && y*scale < templateSizes[i][1])
					{
						avgImg.setColor(x, y, tmp(x*scale, y*scale)/noMatches);
					}
				}
			}
		}
	}
	
	Img result(orgTempSize[0], orgTempSize[1]);
	float scale = (float)(width)/(float)(orgTempSize[0]);
	
	for(int x=0;x<orgTempSize[0];x++)
	{
		for(int y=0; y<orgTempSize[1]; y++)
		{
			result.setColor(x, y, avgImg(x*scale, y*scale));
		}
	}
	
	return result;
	
}

void RepetitionFinder::performPCA(vector<Img *> &templateImg, vector<vector<repetitionGrid> > &matches, int gridIndex, vector<Vec2i> &templateSizes)
{
	cv::Mat Mcomp;
	
	int noImages = templateImg.size();
	int noMatches = 0;
	int noTempMatches = 0;
	int width, height;
	
	//find min width/height
	for(int i=0; i<noImages; i++)
	{
		if(i==0 && templateSizes[i][0] != 0)
		{
			width = templateSizes[i][0];
			height = templateSizes[i][1];
		}
		else if(templateSizes[i][0] < width && templateSizes[i][0] != 0)
		{
			width = templateSizes[i][0];
			height = templateSizes[i][1];
		}
		/*if(i==0 && matches[i].xTrans[0] != 0.0 && matches[i].yTrans[1] != 0.0)
		{
			width = matches[i].xTrans[0]*2;
			height = matches[i].yTrans[1]*2;
		}
		else if(matches[i].xTrans[0] < width && matches[i].xTrans[0] != 0.0 && matches[i].yTrans[1] != 0.0)
		{
			width = matches[i].xTrans[0]*2;
			height = matches[i].yTrans[1]*2;
		}*/
		for(int m=0; m<matches[i][gridIndex].gridCells.size(); m++)
		{
			if(matches[i][gridIndex].gridCells[m].temporary)
				noTempMatches += 1;
			else
				noMatches += 1;
		}
		//noMatches += matches[i].gridCells.size();
	}
	
	//cv::Mat M = cv::Mat(noMatches, width*height, CV_32F);
	cv::Mat MR = cv::Mat(noMatches, width*height, CV_32F);
	cv::Mat MG = cv::Mat(noMatches, width*height, CV_32F);
	cv::Mat MB = cv::Mat(noMatches, width*height, CV_32F);
	cv::Mat Mtemp = cv::Mat(noTempMatches, width*height, CV_32F);
	cv::Mat mean = cv::Mat(1, width*height, CV_32F, 0.0);
	cv::Mat meanR = cv::Mat(1, width*height, CV_32F, 0.0);
	cv::Mat meanG = cv::Mat(1, width*height, CV_32F, 0.0);
	cv::Mat meanB = cv::Mat(1, width*height, CV_32F, 0.0);
	
	float* meanI = mean.ptr<float>(0);
	float* meanIR = meanR.ptr<float>(0);
	float* meanIG = meanG.ptr<float>(0);
	float* meanIB = meanB.ptr<float>(0);
	
	int count = 0;
	int countTemp = 0;
	
	for(int i = 0; i < noImages; i++)
	{
		float scale = (float)(templateSizes[i][0]) / float(width);
		//float scale = matches[i].xTrans[0]*2.0 / (float) width;
		for(int m=0; m<matches[i][gridIndex].gridCells.size(); m++)
		{
			Img tmp(templateSizes[i][0], templateSizes[i][1]);
			//int largerWidth = matches[i].xTrans[0]*2;
			//int largerHeight = matches[i].yTrans[1]*2;
			//Img tmp(largerWidth, largerHeight);
			for(int x=0; x</*largerWidth*/templateSizes[i][0]; x++)
			{
				for(int y=0; y</*largerHeight*/templateSizes[i][1]; y++)
				{
					//if(matches[i].gridCells[m].center[0]-largerWidth/2+x >= 0 && matches[i].gridCells[m].center[0]-largerWidth/2+x < templateImg[i]->width() &&
					  // matches[i].gridCells[m].center[1]-largerHeight/2+y >= 0 && matches[i].gridCells[m].center[1]-largerHeight/2+y < templateImg[i]->height())
					{
						tmp.setColor(x, y, (*templateImg[i])(matches[i][gridIndex].gridCells[m].center[0]-templateSizes[i][0]/2/*largerWidth/2*/+x, 
															 matches[i][gridIndex].gridCells[m].center[1]-templateSizes[i][1]/2/*largerHeight/2*/+y));
					}
					//else
					//{
					//	tmp.setColor(x, y, Color(0.0,0.0,0.0));
					//}
				}
			}
			//tmp.write("/Users/ceylan/Desktop/tmp1.png");
			Img scaledTmp(width, height);
			for(int x=0; x<width; x++)
			{
				for(int y=0; y<height; y++)
				{
					if(x*scale>=0 && x*scale < templateSizes[i][0]/*largerWidth*/ && y*scale>=0 && y*scale < templateSizes[i][1]/*largerHeight*/)
						scaledTmp.setColor(x, y, tmp(x*scale, y*scale));
					else
						scaledTmp.setColor(x, y, Color(0.0,0.0,0.0));
				}
			}
			//scaledTmp.write("/Users/ceylan/Desktop/tmp2.png");
			
			if(!matches[i][gridIndex].gridCells[m].temporary)
			{
				Img tmp;
				tmp.resize(width, height);
				
				//float* Mi = M.ptr<float>(count);
				float* Mir = MR.ptr<float>(count);
				float* Mig = MG.ptr<float>(count);
				float* Mib = MB.ptr<float>(count);
				
				for(int j = 0; j < MR.cols; j++)
				{
					int w = j%width;
					int h = j/width;
					Color c = (scaledTmp)(w, h);
					tmp.setColor(w, h, c);
					
					//Mi[j] = (c[0]+c[1]+c[2])/3.0;
					Mir[j] = c[0];
					Mig[j] = c[1];
					Mib[j] = c[2];
					//meanI[j] += (c[0]+c[1]+c[2])/3.0/(float)noMatches;
					meanIR[j] += c[0]/(float)noMatches;
					meanIG[j] += c[1]/(float)noMatches;
					meanIB[j] += c[2]/(float)noMatches;
				}
				stringstream name;
				name << "/Users/ceylan/Desktop/modes/tmp" << count << ".png";
				string n;
				name >> n;
				tmp.write(n.c_str());
				
				count += 1;
			}
			else
			{
				float* Mi = Mtemp.ptr<float>(countTemp);
				
				for(int j = 0; j < Mtemp.cols; j++)
				{
					int w = j%width;
					int h = j/width;
					Color c = (scaledTmp)(w, h);
					Mi[j] = (c[0]+c[1]+c[2])/3.0;
				}
				countTemp += 1;				
			}
		}
	}
	
	int maxComponents = 20;
	
	//cv::PCA *pca = new cv::PCA(M, mean, CV_PCA_DATA_AS_ROW,  maxComponents);
	cv::PCA *pcaR = new cv::PCA(MR, meanR, CV_PCA_DATA_AS_ROW,  maxComponents);
	cv::PCA *pcaG = new cv::PCA(MG, meanG, CV_PCA_DATA_AS_ROW,  maxComponents);
	cv::PCA *pcaB = new cv::PCA(MB, meanB, CV_PCA_DATA_AS_ROW,  maxComponents);
	
	/*cv::Mat modes = pca->eigenvectors;
	for(int i=0; i<maxComponents; i++)
	{
		float *modesI = modes.ptr<float>(i);
		Img a(width, height);
		for(int j=0; j<modes.cols; j++)
		{
			int w = j%width;
			int h = j/width;
			a.setColor(w,h, Color(modesI[j], modesI[j], modesI[j]));
		}
		stringstream ss;
		ss << "/Users/ceylan/Desktop/modes/" << i << ".png";
		string filename;
		ss >> filename;
		a.write(filename.c_str());
	}*/
	
	count = 0;
	countTemp = 0;
	ofstream fout("/Users/ceylan/Desktop/data.txt", ios::out);
	ofstream fout2("/Users/ceylan/Desktop/data2.txt", ios::out);
	vector<float> projValuesR;
	vector<float> projValuesG;
	vector<float> projValuesB;
	projValuesR.resize(maxComponents, 0.0);
	projValuesG.resize(maxComponents, 0.0);
	projValuesB.resize(maxComponents, 0.0);
	
	for(int i=0; i<noImages; i++)
	{
		//printf("image %d\n", i);
		for(int m=0; m<matches[i][gridIndex].gridCells.size(); m++)
		{
			cv::Mat tmpMatR, tmpMatG, tmpMatB;
			cv::Mat MresR, MresG, MresB;
			if(!matches[i][gridIndex].gridCells[m].temporary)
			{
				tmpMatR = MR.row(count);
				tmpMatG = MG.row(count);
				tmpMatB = MB.row(count);
			}
			//else
			//	tmpMat = Mtemp.row(countTemp);
			
			cv::Mat compR = pcaR->project(tmpMatR);
			cv::Mat compG = pcaG->project(tmpMatG);
			cv::Mat compB = pcaB->project(tmpMatB);
			
			float *compIR = compR.ptr<float>(0);
			float *compIG = compG.ptr<float>(0);
			float *compIB = compB.ptr<float>(0);
			
			if(matches[i][gridIndex].gridCells[m].rowNumber == 0 && matches[i][gridIndex].gridCells[m].colNumber == 0)
				fout << compIR[0] << " " << compIR[1] << endl;
			if(i==0)
			{
				fout2 << compIR[0] << " " << compIR[1] << endl;
				if(matches[i][gridIndex].gridCells[m].colNumber == 2 || matches[i][gridIndex].gridCells[m].colNumber == 1)
				{
					for(int m=0; m<maxComponents; m++)
					{
						projValuesR[m] += compIR[m];
						projValuesG[m] += compIG[m];
						projValuesB[m] += compIB[m];
					}
				}
			}
					
			MresR = pcaR->backProject(compR);
			MresG = pcaG->backProject(compG);
			MresB = pcaB->backProject(compB);
			
			float *MresIR = MresR.ptr<float>(0);
			float *MresIG = MresG.ptr<float>(0);
			float *MresIB = MresB.ptr<float>(0);
			
			Img back;
			back.resize(width, height);
			for(int j=0; j<MresR.cols; j++)
			{
				int w = j%width;
				int h = j/width;
				back.setColor(w,h, Color(MresIR[j], MresIG[j], MresIB[j]));
			}
			stringstream ss;
			ss << "/Users/ceylan/Desktop/modes/back" << count << ".png";
			string filename;
			ss >> filename;
			back.write(filename.c_str());
			
			float scale = (float)(templateSizes[i][0]) / float(width);
			for(int x=0; x<templateSizes[i][0]; x++)
			{
				for(int y=0; y<templateSizes[i][1]; y++)
				{
					{
						int newX = x+matches[i][gridIndex].gridCells[m].center[0]-templateSizes[i][0]/2;
						int newY = y+matches[i][gridIndex].gridCells[m].center[1]-templateSizes[i][1]/2;
						if(newX >=0 && newX<templateImg[i]->width() && newY >=0 && newY<templateImg[i]->height())
						{
							(*templateImg[i]).setColor(x+matches[i][gridIndex].gridCells[m].center[0]-templateSizes[i][0]/2,
													   y+matches[i][gridIndex].gridCells[m].center[1]-templateSizes[i][1]/2,
													   back(x/scale, y/scale));
						}
					}
				}
			}
			
			
			if(!matches[i][gridIndex].gridCells[m].temporary)
				count+=1;
			else
				countTemp+=1;
			
		}
		
		stringstream sAll;
		sAll << "/Users/ceylan/Desktop/modes/all" << i << ".png";
		string filenameAll;
		sAll >> filenameAll;
		(*templateImg[i]).write(filenameAll.c_str());
		
	}
	
	fout.close();
	
	cv::Mat modesR = pcaR->eigenvectors;
	cv::Mat modesG = pcaG->eigenvectors;
	cv::Mat modesB = pcaB->eigenvectors;
	

	Img meanImg(width, height);
	Img interpolated(width, height);
	for(int j=0; j<mean.cols; j++)
	{
		int w = j%width;
		int h = j/width;
		Color c(meanIR[j], meanIG[j], meanIB[j]);
		meanImg.setColor(w,h, c);
		interpolated.setColor(w, h, c);
	}
	for(int m=0; m<maxComponents; m++)
	{
		float *modesIR = modesR.ptr<float>(m);
		float *modesIG = modesG.ptr<float>(m);
		float *modesIB = modesB.ptr<float>(m);
		for(int j=0; j<modesR.cols; j++)
		{
			int w = j%width;
			int h = j/width;
			Color c = interpolated(w, h);
			c[0] = c[0] + modesIR[j]*(projValuesR[m]/2.0);
			c[1] = c[1] + modesIG[j]*(projValuesG[m]/2.0);
			c[2] = c[2] + modesIB[j]*(projValuesB[m]/2.0);
			interpolated.setColor(w,h,c);
		}
	}
	
	meanImg.write("/Users/ceylan/Desktop/mean.png");
	interpolated.write("/Users/ceylan/Desktop/interpolated.png");
}

float RepetitionFinder::computeNCCScoreWithGrayScale(vector<Color> &refPatch, vector<Color> &neighborPatch)
{
	float score = 0.0;
	
	float mean1 = 0.0;
	float mean2 = 0.0;
	float stdDev1 = 0.0;
	float stdDev2 = 0.0;
	float r, n;
	
	int pixelCount = refPatch.size();
	int pixelCount2 = neighborPatch.size();
	if(pixelCount != pixelCount2)
	{
		printf("The number of pixels in two patches don't match.\n");
		return -1.0;
	}
	
	for(int i=0; i<pixelCount; i++)
	{
		r = (refPatch[i][0] + refPatch[i][1] + refPatch[i][2]) / 3.0;
		mean1 += r;
		stdDev1 += r * r;
		
		n = (neighborPatch[i][0] + neighborPatch[i][1] + neighborPatch[i][2]) / 3.0;
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
		r = (refPatch[i][0] + refPatch[i][1] + refPatch[i][2]) / 3.0;
		n = (neighborPatch[i][0] + neighborPatch[i][1] + neighborPatch[i][2]) / 3.0;
		
		score += (r-mean1)*(n-mean2);
		//if(i==(int)(pixelCount/2) && score < (pixelCount/2) * stdDev1 * stdDev2 * 0.4)
		//return -1.0;
	}
	
	score /= (stdDev1*stdDev2*pixelCount);
	//printf("s:%f\n", score);
	return score;
}

float RepetitionFinder::computeNCCScoreWithGrayScale(Vec2f refPixel, Vec2f neighborPixel, Photo &refPhoto, Photo &neighPhoto, float windowSize, int level)
{
	vector<Color> refPatch;
	vector<Color> neighborPatch;
	
	Vec2f xAxis(1.0, 0.0);
	Vec2f yAxis(0.0, 1.0);
	int normalizef = 0;
	
	refPhoto.grabTexture(level, refPixel,xAxis, yAxis, windowSize, refPatch, normalizef);
	neighPhoto.grabTexture(level, neighborPixel, xAxis, yAxis, windowSize, neighborPatch, normalizef);
	
	float pixelCount = float(windowSize * windowSize);
	if(refPatch.size() < (windowSize*windowSize) || neighborPatch.size() < (windowSize*windowSize))
		return 0.0;
	
	return computeNCCScoreWithGrayScale(refPatch, neighborPatch);
}

int RepetitionFinder::findBestMatch(Vec2f refPixel, vector<Vec2f> &candidateMatches, Photo &refPhoto, Photo &neighPhoto, int level)
{
	float windowSize = 11;
	int noMatches = candidateMatches.size();
	int maxScoreIndex = -1;
	float maxScore;
	for(int i=0; i<candidateMatches.size(); i++)
	{
		float score = computeNCCScoreWithGrayScale(refPixel, candidateMatches[i], refPhoto, neighPhoto, windowSize, 0);
		if(maxScoreIndex == -1)
		{
			maxScoreIndex = i;
			maxScore = score;
		}
		else if(score > maxScore)
		{
			maxScore = score;
			maxScoreIndex = i;
		}
	}
	
	if(maxScore < 0.6)
		return -1;
	else
		return maxScoreIndex;
}

void RepetitionFinder::templateMatchingWithSmallerFFTCorrelation(Img* templateImg, Img* sourceImg, vector<Vec2f> &matches, vector<float> &finalScores, float thresh)
{
	int sourceWidth, sourceHeight;
	int templateWidth, templateHeight;
	
	sourceWidth = sourceImg->width();
	sourceHeight = sourceImg->height();
	
	templateWidth = templateImg->width();
	templateHeight = templateImg->height();
	
	int smallWindowSize = 15;//min(templateWidth, templateHeight)/4;
	
	//divide the image into smaller windows
	int noHorWindows = templateWidth/smallWindowSize;
	int noVerWindows = templateHeight/smallWindowSize;
	
	thresh = (0.8)*(0.8)*noHorWindows*noVerWindows;
	
	int tempStartX = (templateWidth - noHorWindows*smallWindowSize)/2 ;
	int tempStartY = (templateHeight - noVerWindows*smallWindowSize)/2 ;
	
	templateWidth = noHorWindows*smallWindowSize;
	templateHeight = noVerWindows*smallWindowSize;
	
	vector<vector<vector<Color> > > templatePatches;
	templatePatches.resize(noHorWindows);
	for(int i=0; i<templatePatches.size(); i++)
		templatePatches[i].resize(noVerWindows);
	
	IplImage *sourceImageCV;
	sourceImageCV = cvCreateImage(cvSize(sourceWidth,sourceHeight),IPL_DEPTH_8U,3);
	for(unsigned int y=0;y<sourceHeight;y++)
	{
		for(unsigned int x=0;x<sourceWidth;x++)
		{	
			CvScalar s;
			s.val[0]= (*sourceImg)(x,y)[2] * 255.0;
			s.val[1]= (*sourceImg)(x,y)[1] * 255.0;
			s.val[2]= (*sourceImg)(x,y)[0] * 255.0;
			
			cvSet2D(sourceImageCV,y,x,s); // set the (i,j) pixel value
		}
	}
	
	float *scores = new float[sourceWidth*sourceHeight];
	memset(scores, 0.0, sizeof(float)*sourceWidth*sourceHeight);
	
	int halfTempWidth = templateWidth/2;
	int halfTempHeight = templateHeight/2;
	
	int halfSmallWidth = smallWindowSize/2;
	int halfSmallHeight = smallWindowSize/2;
	
	//prepare small template patches
//#pragma omp parallel for
	for(int i=0; i<noHorWindows; i++)
	{
		for(int j=0; j<noVerWindows; j++)
		{
			IplImage *templateImageCV, *result;
			
			templateImageCV = cvCreateImage(cvSize(smallWindowSize,smallWindowSize),IPL_DEPTH_8U,3);
			for(unsigned int y=0;y<smallWindowSize;y++)
			{
				for(unsigned int x=0;x<smallWindowSize;x++)
				{	
					CvScalar s;
					s.val[0]= (*templateImg)(tempStartX+i*smallWindowSize+x,tempStartY+j*smallWindowSize+y)[2] * 255.0;
					s.val[1]= (*templateImg)(tempStartX+i*smallWindowSize+x,tempStartY+j*smallWindowSize+y)[1] * 255.0;
					s.val[2]= (*templateImg)(tempStartX+i*smallWindowSize+x,tempStartY+j*smallWindowSize+y)[0] * 255.0;
					
					cvSet2D(templateImageCV,y,x,s); // set the (i,j) pixel value
				}
			}
			
			CvSize resultSize;
			resultSize.width = sourceWidth - smallWindowSize + 1;
			resultSize.height = sourceHeight - smallWindowSize + 1;
			
			result = cvCreateImage( resultSize, IPL_DEPTH_32F, 1 );
			cvMatchTemplate(sourceImageCV, templateImageCV, result, CV_TM_CCOEFF_NORMED);
			
			cvReleaseImage( &templateImageCV );
			
			float* data;
			int step;
			CvSize size;
			
			cvGetRawData( result, ( uchar** )&data, &step, &size );
			
			step /= sizeof( data[0] );
			
			
			for( int y = 0; y < size.height; y++, data += step)
			{
				for( int x = 0; x < size.width; x++ )
				{
					int posX = x + halfSmallWidth;
					int posY = y + halfSmallHeight;
					int scoreX = posX + halfTempWidth - halfSmallWidth - i*smallWindowSize;
					int scoreY = posY + halfTempHeight - halfSmallHeight - j*smallWindowSize;
					if(scoreX>=halfTempWidth && scoreX < sourceWidth-halfTempWidth && scoreY>=halfTempHeight && scoreY<sourceHeight-halfTempHeight)
						scores[scoreY*sourceWidth+scoreX] += (1.0 - data[x])*(1.0-data[x]); // / (float)(noHorWindows*noVerWindows);
				}
			}
			cvReleaseImage( &result );
			
			/*for(int x=-smallWindowSize/2; x<=smallWindowSize/2; x++)
			{
				for(int y=-smallWindowSize/2; y<smallWindowSize/2; y++)
				{
					templatePatches[i][j].push_back((*templateImg)(tempStartX+i*smallWindowSize+x, tempStartY+j*smallWindowSize+y));
				}
			}*/
		}
	}
	
	
	/*#pragma omp parallel for
	for(int x=halfTempWidth; x<sourceWidth-halfTempWidth; x++)
	{
		for(int y=halfTempHeight; y<sourceHeight-halfTempHeight; y++)
		{
			for(int i=0; i<noHorWindows; i++)
			{
				for(int j=0; j<noVerWindows; j++)
				{
					//prepate the source patch
					vector<Color> sourcePatch;
					for(int a=-smallWindowSize/2; a<=smallWindowSize/2; a++)
					{
						for(int b=-smallWindowSize/2; b<smallWindowSize/2; b++)
						{
							sourcePatch.push_back((*sourceImg)(x-halfTempWidth+halfSmallWidth+i*smallWindowSize+a, y-halfTempHeight+halfSmallHeight+j*smallWindowSize+b));
						}
					}
					
					scores[y*sourceWidth+x] += computeNCCScoreWithGrayScale(templatePatches[i][j], sourcePatch);
				}
			}
			
			scores[y*sourceWidth+x] /= (float)(noHorWindows*noVerWindows);
		}
	}*/
	
	Img r(sourceWidth, sourceHeight);
	for(int i=0; i<sourceWidth; i++)
	{
		for(int j=0; j<sourceHeight; j++)
		{
			float s = scores[j*sourceWidth+i]/(4.0*noHorWindows*noVerWindows); 
			r.setColor(i, j, Color(s, s, s));
		}
	}
	r.write("/Users/ceylan/Desktop/result.png");
	
	float *localMaxData = new float[sourceWidth*sourceHeight];
	int halfWindowSize = min(templateWidth, templateHeight) / 2;
	Img scoreImg(sourceWidth, sourceHeight);
	#pragma omp parallel for
	for( int y = 0; y < sourceHeight; y++)
	{
		for( int x = 0; x < sourceWidth; x++ )
		{
			scoreImg.setColor(x, y, Color((scores[y*sourceWidth+x]+1.0)/2.0, (scores[y*sourceWidth+x]+1.0)/2.0, (scores[y*sourceWidth+x]+1.0)/2.0));
			bool localMaxima = true;
			if(x<halfTempWidth || x >= sourceWidth-halfTempWidth || y < halfTempHeight || y >= sourceHeight-halfTempHeight)
			{
				localMaxData[y*sourceWidth+x] = 4.0*noHorWindows*noVerWindows;
				continue;
			}
				
			for(int i=-halfWindowSize; i<=halfWindowSize; i++)
			{
				for(int j=-halfWindowSize; j<=halfWindowSize; j++)
				{
					if(x+i>=halfTempWidth && x+i<sourceWidth-halfTempWidth && y+j>=halfTempHeight && y+j<sourceHeight-halfTempHeight)
					{
						if(scores[(y+j)*sourceWidth + x+i] < scores[y*sourceWidth+x])
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
				localMaxData[y*sourceWidth+x] = scores[y*sourceWidth+x];
				//printf("localMax:%f\n",scores[y*sourceWidth+x]);
			}
			else
			{
				localMaxData[y*sourceWidth+x] = 4.0*noHorWindows*noVerWindows;
			}
		}
	}
	
	scoreImg.write("/Users/ceylan/Desktop/sonuc.png");
	
	for(int y=0; y<sourceHeight; y++)
	{
		for(int x=0; x<sourceWidth; x++)
		{
			if(localMaxData[y*sourceWidth+x] <= thresh)
			{
				//count the number of black pixels
				/*int black = 0;
				for(int w=-halfTempWidth; w<=templateWidth; w++)
				{
					for(int h=-halfTempHeight; h<=templateHeight; h++)
					{
						if((*sourceImg)(x+w, y+h) == Color(0,0,0))
							black++;
					}
				}
				
				if(black > templateWidth*templateHeight/2.0)
					continue;
				*/
				printf("(%d, %d):score:%f\n", x, y, localMaxData[y*sourceWidth+x]);
				matches.push_back(Vec2f(x, y));
				finalScores.push_back(localMaxData[y*sourceWidth+x]);
			}
		}
	}
	
	delete [] localMaxData;
	delete [] scores;
}

void RepetitionFinder::findRepetitionInROI(Img *sourceImg, Img* templateImg, Vec2f upperCorner, Vec2f lowerCorner, vector<Vec2f> &matches, vector<float> &scores, float thresh)
{
	//create source image with search region
	if(upperCorner[0] < 0)
		upperCorner[0] = 0;
	if(lowerCorner[0] > sourceImg->width())
		lowerCorner[0] = sourceImg->width();
	
	if(upperCorner[1] < 0)
		upperCorner[1] = 0;
	else if(lowerCorner[1] > sourceImg->height())
		lowerCorner[1] = sourceImg->height();
	
	int width = lowerCorner[0] - upperCorner[0];
	int height = lowerCorner[1] - upperCorner[1];
	
	if(width < templateImg->width() || height < templateImg->height())
	{
		return;
	}
	
	Img newSourceImg(width, height);
	for(int x=0; x<width; x++)
	{
		for(int y=0; y<height; y++)
		{
			newSourceImg.setColor(x, y, (*sourceImg)(x+upperCorner[0], y+upperCorner[1]));
		}
	}
	
	templateMatchingWithFFTCorrelation(templateImg, &newSourceImg, matches, scores, thresh);
	
	for(int i=0; i<matches.size(); i++)
	{
		matches[i] += upperCorner;
	}
}

void RepetitionFinder::templateMatchingWithFFTCorrelation(Img* templateImg, Img* sourceImg, vector<Vec2f> &matches, vector<float> &scores, float thresh)
{
	int sourceWidth, sourceHeight;
	
	int templateWidth, templateHeight;
	
	IplImage *sourceImageCV, *templateImageCV, *result;
	
	templateWidth = templateImg->width();
	templateHeight = templateImg->height();
	
	sourceWidth = sourceImg->width();
	sourceHeight = sourceImg->height();
	
	if(templateWidth <= 0 || templateHeight <= 0 || sourceWidth <= 0 || sourceHeight <= 0)
		return ;
	
	if(templateWidth > sourceWidth || templateHeight > sourceHeight)
		return ;
	
	templateImageCV = cvCreateImage(cvSize(templateWidth,templateHeight),IPL_DEPTH_8U,3);
	
#pragma omp parallel for
	for(unsigned int y=0;y<templateHeight;y++)
	{
		for(unsigned int x=0;x<templateWidth;x++)
		{	
			CvScalar s;
			s.val[0]= (*templateImg)(x,y)[2] * 255.0;
			s.val[1]= (*templateImg)(x,y)[1] * 255.0;
			s.val[2]= (*templateImg)(x,y)[0] * 255.0;
			
			cvSet2D(templateImageCV,y,x,s); // set the (i,j) pixel value
		}
	}
	
	sourceImageCV = cvCreateImage(cvSize(sourceWidth,sourceHeight),IPL_DEPTH_8U,3);

#pragma omp parallel for
	for(unsigned int y=0;y<sourceHeight;y++)
	{
		for(unsigned int x=0;x<sourceWidth;x++)
		{	
			CvScalar s;
			s.val[0]= (*sourceImg)(x,y)[2] * 255.0;
			s.val[1]= (*sourceImg)(x,y)[1] * 255.0;
			s.val[2]= (*sourceImg)(x,y)[0] * 255.0;
			
			cvSet2D(sourceImageCV,y,x,s); // set the (i,j) pixel value
		}
	}
	
	
	CvSize resultSize;
	resultSize.width = sourceWidth - templateWidth + 1;
	resultSize.height = sourceHeight - templateHeight + 1;
	
	result = cvCreateImage( resultSize, IPL_DEPTH_32F, 1 );
	//IplImage *corrResult = cvCreateImage(resultSize, IPL_DEPTH_32F, 1); 
	cvMatchTemplate(sourceImageCV, templateImageCV, result, CV_TM_CCOEFF_NORMED);
	//cvMatchTemplate(sourceImageCV, templateImageCV, corrResult, CV_TM_CCORR);
	
	cvReleaseImage( &sourceImageCV );
	cvReleaseImage( &templateImageCV );
	
	//find candidate locations
	//extract the raw data for analysis
	float* data;
	float* localMaxData;// *corrData;
	int step;
	CvSize size;
	
	//cvGetRawData( corrResult, (uchar**) &corrData, &step, &size);
	cvGetRawData( result, ( uchar** )&data, &step, &size );
	
	step /= sizeof( data[0] );
	
	localMaxData = new float[size.height*size.width];
	
	int windowSize = min(templateWidth, templateHeight);
	
	int halfWindowSize = windowSize / 2;
	
	printf("halfWindowSize:%d\n", halfWindowSize);
	
#pragma omp parallel for
	for( int y = 0; y < size.height; y++)
	{
		for( int x = 0; x < size.width; x++ )
		{
			bool localMaxima = true;
			if(data[step*y+x] < -1.0 || data[step*y+x] > 1.0)
			{
				localMaxData[step*y+x] = 0.0;
				continue;
			}
			
			
			for(int i=-halfWindowSize; i<=halfWindowSize; i++)
			{
				for(int j=-halfWindowSize; j<=halfWindowSize; j++)
				{
					if(x+i>=0 && x+i<size.width && y+j>=0 && y+j<size.height)
					{
						float *tmp = data + step*(y + j);
						if(tmp[x+i] > data[step*y+x])
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
				localMaxData[step*y+x] = data[step*y+x];
				if(localMaxData[step*y+x] > 0.3)
					printf("localMax at (%d,%d):%f\n", x, y, localMaxData[step*y+x]);
			}
			else
				localMaxData[step*y+x] = 0.0;
		}
	}
	
	for(int y=0; y<size.height; y++, localMaxData += step /*,corrData += step*/)
	{
		for(int x=0; x<size.width; x++)
		{
			if(localMaxData[x] > thresh)
			{
				//count the number of black pixels
				int black = 0;
				for(int w=0; w<templateWidth; w++)
				{
					for(int h=0; h<templateHeight; h++)
					{
						if((*sourceImg)(x+w, y+h) == Color(0,0,0))
							black++;
					}
				}
				
				if(black > templateWidth*templateHeight/2.0)
					continue;
				
				//printf("(%d, %d):score:%f\n", x+templateWidth/2, y+templateHeight/2, localMaxData[x]);
				matches.push_back(Vec2f((x+templateWidth/2), (y+templateHeight/2)));
				scores.push_back(localMaxData[x]);
			}
		}
	}
	
	localMaxData -= step*(size.height);
	//corrData -= step*size.height;
	delete [] localMaxData;
	//delete [] corrData;
	cvReleaseImage( &result);
	//cvReleaseImage( &corrResult);
}

void RepetitionFinder::templateMatchingWithFFTCorrelation(Img* templateImg, Img* sourceImg, Vec2i offset, Vec2i searchSize, vector<Vec2f> &matches, vector<float> &scores, float thresh)
{
	int templateWidth, templateHeight;
	
	IplImage *sourceImageCV, *templateImageCV, *result;
	
	sourceImageCV = cvCreateImage(cvSize(searchSize[0],searchSize[1]),IPL_DEPTH_8U,3);
	for(unsigned int y=0;y<searchSize[1];y++)
	{
		for(unsigned int x=0;x<searchSize[0];x++)
		{	
			CvScalar s;
			s.val[0]= (*sourceImg)(offset[0]+x,offset[1]+y)[2] * 255.0;
			s.val[1]= (*sourceImg)(offset[0]+x,offset[1]+y)[1] * 255.0;
			s.val[2]= (*sourceImg)(offset[0]+x,offset[1]+y)[0] * 255.0;
			
			cvSet2D(sourceImageCV,y,x,s); // set the (i,j) pixel value
		}
	}
	
	cvSaveImage("/Users/ceylan/Desktop/source.png", sourceImageCV);
	templateWidth = templateImg->width();
	templateHeight = templateImg->height();
	templateImageCV = cvCreateImage(cvSize(templateWidth,templateHeight),IPL_DEPTH_8U,3);
	for(unsigned int y=0;y<templateHeight;y++)
	{
		for(unsigned int x=0;x<templateWidth;x++)
		{	
			CvScalar s;
			s.val[0]= (*templateImg)(x,y)[2] * 255.0;
			s.val[1]= (*templateImg)(x,y)[1] * 255.0;
			s.val[2]= (*templateImg)(x,y)[0] * 255.0;
			
			cvSet2D(templateImageCV,y,x,s); // set the (i,j) pixel value
		}
	}
	
	cvSaveImage("/Users/ceylan/Desktop/target.png", templateImageCV);
	
	CvSize resultSize;
	resultSize.width = searchSize[0] - templateWidth + 1;
	resultSize.height = searchSize[1] - templateHeight + 1;
	
	result = cvCreateImage( resultSize, IPL_DEPTH_32F, 1 );
	//IplImage *corrResult = cvCreateImage(resultSize, IPL_DEPTH_32F, 1); 
	cvMatchTemplate(sourceImageCV, templateImageCV, result, CV_TM_CCOEFF_NORMED);
	//cvMatchTemplate(sourceImageCV, templateImageCV, corrResult, CV_TM_CCORR);
	
	cvReleaseImage( &sourceImageCV );
	cvReleaseImage( &templateImageCV );
	
	//find candidate locations
	//extract the raw data for analysis
	float* data;
	float* localMaxData;// *corrData;
	int step;
	CvSize size;
	
	//cvGetRawData( corrResult, (uchar**) &corrData, &step, &size);
	cvGetRawData( result, ( uchar** )&data, &step, &size );
	
	step /= sizeof( data[0] );
	
	localMaxData = new float[size.height*size.width];
	
	int windowSize = min(templateWidth, templateHeight);
	
	int halfWindowSize = windowSize / 2;
	
	printf("halfWindowSize:%d\n", halfWindowSize);
	
	
	for( int y = 0; y < size.height; y++, data += step , localMaxData += step)
	{
		for( int x = 0; x < size.width; x++ )
		{
			bool localMaxima = true;
			if(data[x] < -1.0 || data[x] > 1.0)
			{
				localMaxData[x] = 0.0;
				continue;
			}
			
			
			for(int i=-halfWindowSize; i<=halfWindowSize; i++)
			{
				for(int j=-halfWindowSize; j<=halfWindowSize; j++)
				{
					if(x+i>=0 && x+i<size.width && y+j>=0 && y+j<size.height)
					{
						float *tmp = data + step*j;
						if(tmp[x+i] > data[x])
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
				localMaxData[x] = data[x];
			}
			else
				localMaxData[x] = 0.0;
		}
	}
	
	localMaxData -= step*(size.height);
	for(int y=0; y<size.height; y++, localMaxData += step /*,corrData += step*/)
	{
		for(int x=0; x<size.width; x++)
		{
			if(localMaxData[x] > thresh)
			{
				printf("(%d, %d):score:%f\n", x+templateWidth/2, y+templateHeight/2, localMaxData[x]);
				matches.push_back(Vec2f((x+templateWidth/2+offset[0]), (y+templateHeight/2+offset[1])));
				scores.push_back(localMaxData[x]);
			}
		}
	}
	
	localMaxData -= step*(size.height);
	//corrData -= step*size.height;
	delete [] localMaxData;
	//delete [] corrData;
	cvReleaseImage( &result);
	//cvReleaseImage( &corrResult);
}

void RepetitionFinder::findGridGeneratingVectors(vector<Vec2f> &matchingCenters, vector<float> &scores, float xThresh, float yThresh, 
												 vector<vector<Vec2f> > &finalGroups, vector<Vec2f> &finalTransformations)
{
	int noCenters = matchingCenters.size();
	vector<Vec2f> tmpMatches;
	tmpMatches.assign(matchingCenters.begin(), matchingCenters.end());
	
	if(noCenters == 0)
	{
		printf("No matches found.\n");
		return ; 
	}
	
	//cluster centers according to y
	vector<int> rowIndices;
	vector<float> yCoords;
	sort(tmpMatches.begin(), tmpMatches.end(), lessWrtYCoord);
	rowIndices.push_back(0);
	
	for(int i=1; i<noCenters; i++)
	{
		if(abs(tmpMatches[i][1]-tmpMatches[i-1][1]) > yThresh)
		{
			float avgY = 0;
			for(int index = rowIndices[rowIndices.size()-1]; index < i; index++)
				avgY += tmpMatches[index][1];
			avgY /= (i-rowIndices[rowIndices.size()-1]);
			yCoords.push_back(avgY);
			rowIndices.push_back(i);
		}
	}
	float avgY = 0;
	for(int index = rowIndices[rowIndices.size()-1]; index < noCenters; index++)
		avgY += tmpMatches[index][1];
	avgY /= (noCenters-rowIndices[rowIndices.size()-1]);
	yCoords.push_back(avgY);
	
	float thresh = 50.0;
	
	vector<gridTransformation> transformationCandidatesY;
	vector<vector<float> > transformations;
	
	for(int i=0; i<yCoords.size(); i++)
	{
		int j= i+1;
		//for(int j=i+1; j<yCoords.size(); j++)
		if(j < yCoords.size())
		{
			float diff = yCoords[j]-yCoords[i];
			int newTransformation = -1;
			for(int t=0; t<transformations.size(); t++)
			{
				if(abs(diff - transformations[t][0]) < thresh)
				{
					transformations[t].push_back(diff);
					float sum = 0;
					for(int c=1; c<transformations[t].size(); c++)
					{
						sum += transformations[t][c];
					}
					transformations[t][0] = sum / (float)(transformations[t].size()-1);
					newTransformation = t;
					break;
				}
				else if(diff < transformations[t][0])
				{
					//if we have T, check for T/2, T/3...
					int factor = floor(transformations[t][0] / diff + 0.5);
					if(abs(diff*factor - transformations[t][0]) < thresh)
					{
						float sum = 0;
						for(int c=1; c<transformations[t].size(); c++)
						{
							transformations[t][c] = transformations[t][c] / factor;
							sum += transformations[t][c];
						}
						transformations[t].push_back(diff);
						sum += diff;
						transformations[t][0] = sum / (float)(transformations[t].size()-1);
						newTransformation = t;
						break;							
					}
				}
				else if(diff > transformations[t][0])
				{
					//if we have T, check for 2T, 3T...
					int factor = floor(diff / transformations[t][0] + 0.5);
					if(abs(diff - factor*transformations[t][0]) < thresh)
					{
						diff /= factor;
						transformations[t].push_back(diff);
						float sum = 0;
						for(int c=1; c<transformations[t].size(); c++)
						{
							sum += transformations[t][c];
						}
						transformations[t][0] = sum / (float)(transformations[t].size()-1);
						newTransformation = t;
						break;
					}
				}				
			}
			if(newTransformation == -1)
			{
				newTransformation = transformations.size();
				vector<float> values;
				values.push_back(diff);
				values.push_back(diff);
				transformations.push_back(values);
				
			}
			gridTransformation gt;
			gt.startIndex = i;
			gt.endIndex = j;
			gt.transformation = Vec2f(0.0, transformations[newTransformation][0]);
			gt.transformationId = newTransformation;
			transformationCandidatesY.push_back(gt);
		}
	}
	
	vector<vector<int> > chainsY;
	for(int i=0; i<transformationCandidatesY.size(); i++)
	{
		printf("chain:%d\n", chainsY.size());
		vector<int> chain;
		chain.push_back(i);
		bool stop = false;
		int currentStart = transformationCandidatesY[i].endIndex;
		int currentTr = transformationCandidatesY[i].transformationId;
		printf("%d to %d,", transformationCandidatesY[i].startIndex, transformationCandidatesY[i].endIndex);
		while(!stop)
		{
			stop = true;
			for(int j=0; j<transformationCandidatesY.size(); j++)
			{
				if(transformationCandidatesY[j].startIndex == currentStart && transformationCandidatesY[j].transformationId == currentTr)
				{
					stop = false;
					chain.push_back(j);
					currentStart = transformationCandidatesY[j].endIndex;
					printf("%d to %d,", transformationCandidatesY[j].startIndex, transformationCandidatesY[j].endIndex);
				}
			}
		}
		chainsY.push_back(chain);
		printf("\n");
	}
	
	vector<vector<Vec2f> > groups;
	vector<Vec2f> yTransformations;
	vector<bool> selected;
	selected.resize(rowIndices.size(), false);
	int noSelected = 0;
	while(noSelected < tmpMatches.size())
	{
		if(chainsY.size() == 0)
			break;
		
		//find longest chain
		int longestChainIndex = -1;
		int longestChainLength = 0;
		for(int i=0; i<chainsY.size(); i++)
		{
			if(i==0)
			{
				longestChainIndex = i;
				longestChainLength = chainsY[i].size();
			}
			else if(chainsY[i].size() > longestChainLength)
			{
				longestChainIndex = i;
				longestChainLength = chainsY[i].size();
			}
		}
	 
		vector<Vec2f> group;
		int lastRow;
		for(int e=0; e<chainsY[longestChainIndex].size(); e++)
		{
			int startIndex = transformationCandidatesY[chainsY[longestChainIndex][e]].startIndex;
			int endIndex = transformationCandidatesY[chainsY[longestChainIndex][e]].endIndex;
			int rowStart = rowIndices[startIndex];
			int rowEnd = rowIndices[endIndex];
			selected[startIndex] = true;
			if(e == chainsY[longestChainIndex].size()-1)
				lastRow = endIndex;
			for(int c=rowStart; c<rowEnd; c++)
			{
				group.push_back(tmpMatches[c]);
				noSelected += 1;
			}
		}
		selected[lastRow] = true;
		int nextRow = lastRow + 1;
		if(nextRow == yCoords.size())
		{
			int rowStart = rowIndices[lastRow];
			for(int c=rowStart; c<tmpMatches.size(); c++)
			{
				group.push_back(tmpMatches[c]);
				noSelected += 1;
			}
		}
		else
		{
			int rowStart = rowIndices[lastRow];
			int rowEnd = rowIndices[nextRow];
			for(int c=rowStart; c<rowEnd; c++)
			{
				group.push_back(tmpMatches[c]);
				noSelected += 1;
			}
		}
		
		groups.push_back(group);
		yTransformations.push_back(Vec2f(0.0,transformations[transformationCandidatesY[chainsY[longestChainIndex][0]].transformationId][0]));
		
		vector<vector<int> > tmpChainY;
		for(int i=0; i<chainsY.size(); i++)
		{
			bool add = true;
			for(int e=0; e<chainsY[longestChainIndex].size(); e++)
			{
				int startIndex = transformationCandidatesY[chainsY[longestChainIndex][e]].startIndex;
				int endIndex = transformationCandidatesY[chainsY[longestChainIndex][e]].endIndex;
				
				for(int j=0; j<chainsY[i].size(); j++)
				{
					if(startIndex == transformationCandidatesY[chainsY[i][j]].startIndex || startIndex == transformationCandidatesY[chainsY[i][j]].endIndex ||
					   endIndex == transformationCandidatesY[chainsY[i][j]].startIndex || endIndex == transformationCandidatesY[chainsY[i][j]].endIndex)
					{
						add = false;
						break;
					}
				}
				if(!add)
					break;
			}
			if(add)
			{
				vector<int> chain;
				for(int e=0; e<chainsY[i].size(); e++)
					chain.push_back(chainsY[i][e]);
				tmpChainY.push_back(chain);
			}
		}
		
		chainsY.clear();
		chainsY = tmpChainY;
	}
	for(int i=0; i<rowIndices.size(); i++)
	{
		if(!selected[i])
		{
			//add a new group
			vector<Vec2f> group;
			int rowStart = rowIndices[i];
			int rowEnd;
			if(i+1 == yCoords.size())
				rowEnd = tmpMatches.size();
			else
				rowEnd = rowIndices[i+1];
			for(int c=rowStart; c<rowEnd; c++)
				group.push_back(tmpMatches[c]);
			groups.push_back(group);
			yTransformations.push_back(Vec2f(0.0, 0.0));
		}
	}
	
	//cluster centers according to x
	for(int g=0; g<groups.size(); g++)
	{
		vector<int> colIndices;
		vector<float> xCoords;
		//sort(tmpMatches.begin(), tmpMatches.end(), lessWrtXCoord);
		sort(groups[g].begin(), groups[g].end(), lessWrtXCoord);
		colIndices.push_back(0);
		noCenters = groups[g].size();
		for(int i=1; i<noCenters; i++)
		{
			//if(abs(tmpMatches[i][0]-tmpMatches[i-1][0]) > xThresh)
			if(abs(groups[g][i][0]-groups[g][i-1][0]) > xThresh)
			{
				float avgX = 0;
				for(int index = colIndices[colIndices.size()-1]; index < i; index++)
					//avgX += tmpMatches[index][0];
					avgX += groups[g][index][0];
				avgX /= (i-colIndices[colIndices.size()-1]);
				xCoords.push_back(avgX);
				colIndices.push_back(i);
			}
		}
		float avgX = 0;
		for(int index = colIndices[colIndices.size()-1]; index < noCenters; index++)
			avgX += groups[g][index][0];
		avgX /= (noCenters-colIndices[colIndices.size()-1]);
		xCoords.push_back(avgX);
	
		vector<gridTransformation> transformationCandidatesX;
		for(int i=0; i<transformations.size(); i++)
			transformations[i].clear();
		transformations.clear();
	
		for(int i=0; i<xCoords.size(); i++)
		{
			int j= i+1;
			//for(int j=i+1; j<xCoords.size(); j++)
			if(j < xCoords.size())
			{
				float diff = xCoords[j]-xCoords[i];
				int newTransformation = -1;
				for(int t=0; t<transformations.size(); t++)
				{
					if(abs(diff - transformations[t][0]) < thresh)
					{
						transformations[t].push_back(diff);
						float sum = 0;
						for(int c=1; c<transformations[t].size(); c++)
						{
							sum += transformations[t][c];
						}
						transformations[t][0] = sum / (float)(transformations[t].size()-1);
						newTransformation = t;
						break;
					}
					else if(diff < transformations[t][0])
					{
						//if we have T, check for T/2, T/3...
						int factor = floor(transformations[t][0] / diff + 0.5);
						if(abs(diff*factor - transformations[t][0]) < thresh)
						{
							float sum = 0;
							for(int c=1; c<transformations[t].size(); c++)
							{
								transformations[t][c] = transformations[t][c] / factor;
								sum += transformations[t][c];
							}
							transformations[t].push_back(diff);
							sum += diff;
							transformations[t][0] = sum / (float)(transformations[t].size()-1);
							newTransformation = t;
							break;							
						}
					}
					else if(diff > transformations[t][0])
					{
						//if we have T, check for 2T, 3T...
						int factor = floor(diff / transformations[t][0] + 0.5);
						if(abs(diff - transformations[t][0]*factor) < thresh)					
						{
							printf("transformations[t]:%f\n", transformations[t][0]);
							diff /= factor;
							transformations[t].push_back(diff);
							float sum = 0;
							for(int c=1; c<transformations[t].size(); c++)
							{
								sum += transformations[t][c];
							}
							transformations[t][0] = sum / (float)(transformations[t].size()-1);
							newTransformation = t;
							break;
						}
					}				
				}
				if(newTransformation == -1)
				{
					printf("transformations:%f\n", diff);
					newTransformation = transformations.size();
					vector<float> values;
					values.push_back(diff);
					values.push_back(diff);
					transformations.push_back(values);
				}
			
				gridTransformation gt;
				gt.startIndex = i;
				gt.endIndex = j;
				gt.transformation = Vec2f(transformations[newTransformation][0], 0.0);
				gt.transformationId = newTransformation;
				transformationCandidatesX.push_back(gt);
			}
		}
	
		vector<vector<int> > chainsX;
		for(int i=0; i<transformationCandidatesX.size(); i++)
		{
			printf("chain:%d\n", chainsX.size());
			vector<int> chain;
			chain.push_back(i);
			bool stop = false;
			int currentStart = transformationCandidatesX[i].endIndex;
			int currentTr = transformationCandidatesX[i].transformationId;
			printf("%d to %d,", transformationCandidatesX[i].startIndex, transformationCandidatesX[i].endIndex);
			while(!stop)
			{
				stop = true;
				for(int j=0; j<transformationCandidatesX.size(); j++)
				{
					if(transformationCandidatesX[j].startIndex == currentStart && transformationCandidatesX[j].transformationId == currentTr)
					{
						stop = false;
						chain.push_back(j);
						currentStart = transformationCandidatesX[j].endIndex;
						printf("%d to %d,", transformationCandidatesX[j].startIndex, transformationCandidatesX[j].endIndex);
					}
				}
			}
			chainsX.push_back(chain);
			printf("\n");
		}
	
		int noSelected = 0;
		vector<bool> selected;
		selected.resize(xCoords.size(), false);
		while(noSelected < groups[g].size())
		{
			if(chainsX.size() == 0)
				break;
			
			//find longest chain
			int longestChainIndex = -1;
			int longestChainLength = 0;
			for(int i=0; i<chainsX.size(); i++)
			{
				if(i==0)
				{
					longestChainIndex = i;
					longestChainLength = chainsX[i].size();
				}
				else if(chainsX[i].size() > longestChainLength)
				{
					longestChainIndex = i;
					longestChainLength = chainsX[i].size();
				}
			}
		
			vector<Vec2f> group;
			int lastCol;
			for(int e=0; e<chainsX[longestChainIndex].size(); e++)
			{
				int startIndex = transformationCandidatesX[chainsX[longestChainIndex][e]].startIndex;
				int endIndex = transformationCandidatesX[chainsX[longestChainIndex][e]].endIndex;
				selected[startIndex] = true;
				if(e==chainsX[longestChainIndex].size()-1)
					lastCol = endIndex;
				int colStart = colIndices[startIndex];
				int colEnd = colIndices[endIndex];
				for(int c=colStart; c<colEnd; c++)
				{
					group.push_back(groups[g][c]);
					noSelected += 1;
				}
			}
			selected[lastCol] = true;
			int nextCol = lastCol + 1;
			if(nextCol == xCoords.size())
			{
				int colStart = colIndices[lastCol];
				for(int c=colStart; c<groups[g].size(); c++)
				{
					group.push_back(groups[g][c]);
					noSelected += 1;
				}
			}
			else
			{
				int colStart = colIndices[lastCol];
				int colEnd = colIndices[nextCol];
				for(int c=colStart; c<colEnd; c++)
				{
					group.push_back(groups[g][c]);
					noSelected += 1;
				}
			}
			finalGroups.push_back(group);
			finalTransformations.push_back(Vec2f(transformations[transformationCandidatesX[chainsX[longestChainIndex][0]].transformationId][0], yTransformations[g][1]));
		
			vector<vector<int> > tmpChainX;
			for(int i=0; i<chainsX.size(); i++)
			{
				bool add = true;
				for(int e=0; e<chainsX[longestChainIndex].size(); e++)
				{
					int startIndex = transformationCandidatesX[chainsX[longestChainIndex][e]].startIndex;
					int endIndex = transformationCandidatesX[chainsX[longestChainIndex][e]].endIndex;
					
					for(int j=0; j<chainsX[i].size(); j++)
					{
						if(startIndex == transformationCandidatesX[chainsX[i][j]].startIndex || startIndex == transformationCandidatesX[chainsX[i][j]].endIndex ||
						   endIndex == transformationCandidatesX[chainsX[i][j]].startIndex || endIndex == transformationCandidatesX[chainsX[i][j]].endIndex)
						{
							add = false;
							break;
						}
					}
					if(!add)
						break;
				}
				if(add)
				{
					vector<int> chain;
					for(int e=0; e<chainsX[i].size(); e++)
						chain.push_back(chainsX[i][e]);
					tmpChainX.push_back(chain);
				}
			}
		
			chainsX.clear();
			chainsX = tmpChainX;
		}
		for(int i=0; i<colIndices.size(); i++)
		{
			if(!selected[i])
			{
				//add a new group
				vector<Vec2f> group;
				int colStart = colIndices[i];
				int colEnd;
				if(i+1 == xCoords.size())
					colEnd = groups[g].size();
				else
					colEnd = colIndices[i+1];
				for(int c=colStart; c<colEnd; c++)
					group.push_back(groups[g][c]);
				finalGroups.push_back(group);
				finalTransformations.push_back(Vec2f(0.0, yTransformations[g][1]));
			}
			
		}
	}
	
	for(int i=0; i<finalGroups.size(); i++)
	{
		printf("group:%d with transformation:%f,%f\n",i, finalTransformations[i][0], finalTransformations[i][1]);
		//for(int j=0; j<finalGroups[i].size(); j++)
		//	printf("%d, " ,finalGroups[i][j]);
		//printf("\n");
	}
}

void RepetitionFinder::segmentIntoRowsAndColumns(vector<Vec2f> &matchingCenters, vector<float> &scores, repetitionGrid &finalGrid, Vec2f gridVectorY, Vec2f gridVectorX, int w, int h)
{
	int noCenters = matchingCenters.size();
	vector<Vec2f> tmpMatches;
	
	if(noCenters == 0)
	{
		printf("No matches found.\n");
		return ;
	}
	
	//remove outliers wrt to y axis
	/*float thresh = gridVectorY[1] / 5.0;
	sort(matchingCenters.begin(), matchingCenters.end(), lessWrtYCoord);
	
	for(int i=0; i<noCenters; i++)
	{
		bool outlier = true;
		for(int j=i+1; j<noCenters; j++)
		{
			if(matchingCenters[j][1]-matchingCenters[i][1] < thresh)
				continue;
			
			float factor = floor((matchingCenters[j][1]-matchingCenters[i][1])/gridVectorY[1] + 0.5);
			if(abs(matchingCenters[j][1]-matchingCenters[i][1]-gridVectorY[1]*factor) < thresh)
			{
				outlier = false;
				break;
			}
		}
		if(outlier)
		{
			for(int j=i-1; j>=0; j--)
			{
				if(matchingCenters[i][1]-matchingCenters[j][1] < thresh)
					continue;
				
				float factor = floor((matchingCenters[i][1]-matchingCenters[j][1])/gridVectorY[1] + 0.5);
				if(abs(matchingCenters[i][1]-matchingCenters[j][1]-gridVectorY[1]*factor) < thresh)
				{
					outlier = false;
					break;
				}
				
			}
		}
		if(!outlier)
		{
			tmpMatches.push_back(matchingCenters[i]);
		}
	}
	
	//remove outliers wrt to x axis
	thresh = gridVectorX[0] / 5.0;
	matchingCenters.clear();
	matchingCenters.assign(tmpMatches.begin(), tmpMatches.end());
	sort(matchingCenters.begin(), matchingCenters.end(), lessWrtXCoord);
	tmpMatches.clear();
	
	for(int i=0; i<noCenters; i++)
	{
		bool outlier = true;
		for(int j=i+1; j<noCenters; j++)
		{
			if(matchingCenters[j][0]-matchingCenters[i][0] < thresh)
				continue;
			
			float factor = floor((matchingCenters[j][0]-matchingCenters[i][0])/gridVectorX[0] + 0.5);
			if(abs(matchingCenters[j][0]-matchingCenters[i][0]-gridVectorX[0]*factor) < thresh)
			{
				outlier = false;
				break;
			}
		}
		if(outlier)
		{
			for(int j=i-1; j>=0; j--)
			{
				if(matchingCenters[i][0]-matchingCenters[j][0] < thresh)
					continue;
				
				float factor = floor((matchingCenters[i][0]-matchingCenters[j][0])/gridVectorX[0] + 0.5);
				if(abs(matchingCenters[i][0]-matchingCenters[j][0]-gridVectorX[0]*factor) < thresh)
				{
					outlier = false;
					break;
				}
				
			}
		}
		if(!outlier)
		{
			tmpMatches.push_back(matchingCenters[i]);
		}
	}
	
	matchingCenters.clear();
	matchingCenters.assign(tmpMatches.begin(), tmpMatches.end());
	*/
	
	noCenters = matchingCenters.size();
	vector<repetitionCell> cells;
	for(int i=0; i<matchingCenters.size(); i++)
	{
		cells.push_back(repetitionCell(matchingCenters[i]));
		cells[i].descriptor = scores[i];
	}
	
	//cluster centers according to y
	vector<int> rowIndices;
	
	sort(cells.begin(), cells.end(), lessWrtCellYCoord);
	rowIndices.push_back(0);
	float yThresh;
	if(gridVectorY[1] == 0.0)
		yThresh = h;
	else
		yThresh = gridVectorY[1] / 2.0;
	
	int rowCount = 0;
	for(int i=1; i<noCenters; i++)
	{
		if(abs(cells[i].center[1]-cells[i-1].center[1]) > yThresh)
		{
			rowIndices.push_back(i);
			for(int j=rowIndices[rowIndices.size()-2]; j<rowIndices[rowIndices.size()-1]; j++)
			{
				cells[j].rowNumber = rowCount;
			}
			rowCount += floor((cells[rowIndices[rowIndices.size()-1]].center[1] - cells[rowIndices[rowIndices.size()-2]].center[1]) / gridVectorY[1] + 0.5);
		}
	}
	
	for(int j=rowIndices[rowIndices.size()-1]; j<noCenters; j++)
	{
		cells[j].rowNumber = rowCount;
	}
	
	finalGrid.numRows = rowCount+1;
	
	//cluster centers according to x
	vector<int> colIndices;
	
	sort(cells.begin(), cells.end(), lessWrtCellXCoord);
	colIndices.push_back(0);
	float xThresh;
	if(gridVectorX[0] == 0.0)
		xThresh = w;
	else
		xThresh = gridVectorX[0] / 2.0;
	
	int colCount = 0;
	for(int i=1; i<noCenters; i++)
	{
		if(abs(cells[i].center[0]-cells[i-1].center[0]) > xThresh)
		{
			colIndices.push_back(i);
			for(int j=colIndices[colIndices.size()-2]; j<colIndices[colIndices.size()-1]; j++)
			{
				cells[j].colNumber = colCount;
			}
			colCount += floor((cells[colIndices[colIndices.size()-1]].center[0] - cells[colIndices[colIndices.size()-2]].center[0]) / gridVectorX[0] + 0.5);
		}
	}
	
	for(int j=colIndices[colIndices.size()-1]; j<noCenters; j++)
	{
		cells[j].colNumber = colCount;
	}
	
	finalGrid.numColumns = colCount+1;
	
	sort(cells.begin(), cells.end(), lessWrtCellYCoord);
	
	for(int i=0; i<rowIndices.size(); i++)
	{
		vector<repetitionCell> xCoords;
		if(i==rowIndices.size()-1)
		{
			for(int j=rowIndices[i]; j<noCenters; j++)
			{
				xCoords.push_back(cells[j]);
			}
		}
		else
		{
			for(int j=rowIndices[i]; j<rowIndices[i+1]; j++)
			{
				xCoords.push_back(cells[j]);
			}
		}
		
		sort(xCoords.begin(), xCoords.end(), lessWrtCellXCoord);
		for(int k=0; k<xCoords.size(); k++)
		{
			finalGrid.gridCells.push_back(xCoords[k]);
			printf("row:%d col:%d score:%f\n\n", finalGrid.gridCells[finalGrid.gridCells.size()-1].rowNumber, finalGrid.gridCells[finalGrid.gridCells.size()-1].colNumber, finalGrid.gridCells[finalGrid.gridCells.size()-1].descriptor);
		}
	}
}

void RepetitionFinder::fillMissingGridCells(repetitionGrid &grid, Img *templateImg, Img *sourceImg, float thresh)
{
	int numRows = grid.numRows;
	int numCols = grid.numColumns;
	
	vector<vector<int> > descriptors;
	descriptors.resize(numRows);
	for(int i=0; i<numRows; i++)
	{
		descriptors[i].resize(numCols, -1);
	}
	for(int i=0; i<grid.gridCells.size(); i++)
	{
		descriptors[grid.gridCells[i].rowNumber][grid.gridCells[i].colNumber] = i;
	}
	
	for(int i=0; i<numRows; i++)
	{
		for(int j=0; j<numCols; j++)
		{
			if(descriptors[i][j] == -1)
			{
				//find closest detection
				int closestIndex = -1;
				float dist;
				for(int x=0; x<numRows; x++)
				{
					for(int y=0; y<numCols; y++)
					{
						if(descriptors[x][y] != -1)
						{
							if(closestIndex == -1)
							{
								closestIndex = x*numCols + y;
								dist = (x-i)*(x-i) + (y-j)*(y-j);
							}
							else if((x-i)*(x-i) + (y-j)*(y-j) < dist)
							{
								dist = (x-i)*(x-i) + (y-j)*(y-j);
								closestIndex = x*numCols + y;
							}
						}
					}
				}
				//generate candidate position
				int r = closestIndex / numCols;
				int c = closestIndex % numCols;
				Vec2f pos = grid.gridCells[descriptors[r][c]].center + grid.xTrans*(j-c) + grid.yTrans*(i-r);
				Vec2f upperCorner(pos[0] - templateImg->width(), pos[1] - templateImg->height());
				Vec2f lowerCorner(pos[0] + templateImg->width(), pos[1] + templateImg->height());
				vector<Vec2f> matches;
				vector<float> scores;
				findRepetitionInROI(sourceImg, templateImg, upperCorner, lowerCorner, matches, scores, thresh);
				if(matches.size() > 1)
				{
					printf("More than one match found!!!\n");
					closestIndex = -1;
					dist = 0;
					for(int m=0; m<matches.size(); m++)
					{
						if(closestIndex == -1)
						{
							closestIndex = m;
							dist = (pos - matches[m]).length();
						}
						else if((pos - matches[m]).length() < dist)
						{
							closestIndex = m;
							dist = (pos - matches[m]).length();
						}
						printf("Adding missing cell at row:%d, col:%d\n", i, j);
						repetitionCell cell(matches[closestIndex]);
						cell.size = Vec2i(templateImg->width(), templateImg->height());
						cell.rowNumber = i;
						cell.colNumber = j;
						cell.descriptor = scores[closestIndex];
						grid.gridCells.push_back(cell);
					}
				}
				else if(matches.size() == 0)
				{
					printf("No match found!!!\n");
				}
				else
				{
					printf("Adding missing cell at row:%d, col:%d\n", i, j);
					repetitionCell cell(matches[0]);
					cell.size = Vec2i(templateImg->width(), templateImg->height());
					cell.rowNumber = i;
					cell.colNumber = j;
					cell.descriptor = scores[0];
					grid.gridCells.push_back(cell);
				}
			}
		}
	}
}

void RepetitionFinder::fillMissingGridCellsTemporarily(repetitionGrid &grid, Vec2i size)
{
	int numRows = grid.numRows;
	int numCols = grid.numColumns;
	
    printf("numRows:%d, numCols:%d\n", grid.numRows, grid.numColumns);
    
	vector<vector<int> > descriptors;
	descriptors.resize(numRows);
	for(int i=0; i<numRows; i++)
	{
		descriptors[i].resize(numCols, -1);
	}
	for(int i=0; i<grid.gridCells.size(); i++)
	{
		descriptors[grid.gridCells[i].rowNumber][grid.gridCells[i].colNumber] = i;
	}
	
	for(int i=0; i<numRows; i++)
	{
		for(int j=0; j<numCols; j++)
		{
			if(descriptors[i][j] == -1)
			{
				//find closest detection
				int closestIndex = -1;
				float dist;
				for(int x=0; x<numRows; x++)
				{
					for(int y=0; y<numCols; y++)
					{
						if(descriptors[x][y] != -1)
						{
							if(closestIndex == -1)
							{
								closestIndex = x*numCols + y;
								dist = (x-i)*(x-i) + (y-j)*(y-j);
							}
							else if((x-i)*(x-i) + (y-j)*(y-j) < dist)
							{
								dist = (x-i)*(x-i) + (y-j)*(y-j);
								closestIndex = x*numCols + y;
							}
						}
					}
				}
				//generate candidate position
				int r = closestIndex / numCols;
				int c = closestIndex % numCols;
				Vec2f pos = grid.gridCells[descriptors[r][c]].center + grid.xTrans*(j-c) + grid.yTrans*(i-r);
				
				repetitionCell cell(pos);
				cell.size = size;
				cell.rowNumber = i;
				cell.colNumber = j;
				cell.descriptor = -1.0;
				cell.temporary = true;
				grid.gridCells.push_back(cell);
				printf("add temp cell at row:%d, col:%d\n", i, j);
			}
		}
	}
}

int RepetitionFinder::translateTransformationToGrid(repetitionGrid &refGrid, repetitionGrid &neighborGrid, float scale, Vec2f translation, 
													float repThresh, int width, int height, int &colShift, int &rowShift)
{
	float dist = repThresh;
	int closestIndex, refIndex;
	closestIndex = -1; refIndex = -1;
	
	int tmpIndex;
	float tmpDist;
	
	int nTmpIndex;
	float nTmpDist;
	
	//printf("**************%d inliers***********\n", inliers.size());
	
	vector<Vec2i> alignments;
	vector<int> alignmentVotes;
	
	int minRow = -1; int minCol = -1;
	tmpIndex = -1;
	for(int c=0; c<neighborGrid.gridCells.size(); c++)
	{
		if(tmpIndex == -1)
		{
			minRow = neighborGrid.gridCells[c].rowNumber;
			minCol = neighborGrid.gridCells[c].colNumber;
			tmpIndex = c;
		}
		else if(neighborGrid.gridCells[c].rowNumber <= minRow && neighborGrid.gridCells[c].colNumber <= minCol)
		{
			minRow = neighborGrid.gridCells[c].rowNumber;
			minCol = neighborGrid.gridCells[c].colNumber;
			tmpIndex = c;
		}
	}
	
	Vec2f featurePos = neighborGrid.gridCells[tmpIndex].center;
	featurePos = featurePos*scale + translation;
		
	nTmpIndex = -1;
		
	//find closest grid center
	for(int c=0; c<refGrid.gridCells.size(); c++)
	{
		if(nTmpIndex == -1)
		{
			nTmpIndex = c;
			nTmpDist = (refGrid.gridCells[c].center-featurePos).length();
		}
		else if((refGrid.gridCells[c].center-featurePos).length() < nTmpDist)
		{
			nTmpIndex = c;
			nTmpDist = (refGrid.gridCells[c].center-featurePos).length();
		}
	}
		
	if(nTmpDist < repThresh)
	{
		rowShift = refGrid.gridCells[nTmpIndex].rowNumber - neighborGrid.gridCells[tmpIndex].rowNumber;
		colShift = refGrid.gridCells[nTmpIndex].colNumber - neighborGrid.gridCells[tmpIndex].colNumber;
		return 1;
	}
	else if(nTmpDist < width && nTmpDist < height)
	{
		//the corresponding cell could be missing
		Vec2f pos = refGrid.gridCells[nTmpIndex].center;
		int noRowsUp = 0;
		int noRowsDown = 0;
		if(refGrid.yTrans[1] != 0.0)
		{
			noRowsUp = floor((float)nTmpDist / (float)refGrid.yTrans[1] + 0.5);
			noRowsDown = floor((float)nTmpDist / (float)refGrid.yTrans[1] + 0.5);
		}
		int noColsLeft = 0;
		int noColsRight = 0;
		if(refGrid.xTrans[0] != 0.0)
		{
			noColsLeft = floor((float)nTmpDist / (float)refGrid.xTrans[0] + 0.5);
			noColsRight = floor((float)nTmpDist / (float)refGrid.xTrans[0] + 0.5);
		}
		int newMinRowIndex = refGrid.gridCells[nTmpIndex].rowNumber;
		int newMinColIndex = refGrid.gridCells[nTmpIndex].colNumber;
		float newMinDistance = nTmpDist;
		for(int y=-noRowsUp; y<=noRowsDown; y++)
		{
			for(int x=-noColsLeft; x<=noColsRight; x++)
			{
				Vec2f newPos = pos + refGrid.yTrans*y + refGrid.xTrans*x;
				if((newPos-featurePos).length() < newMinDistance)
				{
					newMinDistance = (newPos-featurePos).length();
					newMinRowIndex = refGrid.gridCells[nTmpIndex].rowNumber + y;
					newMinColIndex = refGrid.gridCells[nTmpIndex].colNumber + x;
				}
			}
		}
		if(newMinDistance < repThresh)
		{
			rowShift = newMinRowIndex - neighborGrid.gridCells[tmpIndex].rowNumber;
			colShift = newMinColIndex - neighborGrid.gridCells[tmpIndex].colNumber;
			return 1;
		}
		else
			return 0;
	}
}


int RepetitionFinder::translateTransformationToGrid(repetitionGrid &refGrid, repetitionGrid &neighborGrid, vector<Vec2f> &aPts, vector<Vec2f> &bPts, 
													 vector<int> &inliers, float scale, Vec2f translation, float repThresh, int width, int height, int &colShift, int &rowShift)
{
	float dist = repThresh;
	int closestIndex, refIndex;
	closestIndex = -1; refIndex = -1;
	
	int tmpIndex;
	float tmpDist;
	
	int nTmpIndex;
	float nTmpDist;
	
	//printf("**************%d inliers***********\n", inliers.size());
	
	vector<Vec2i> alignments;
	vector<int> alignmentVotes;
	
	int minRow = -1; int minCol = -1;
	tmpIndex = -1;
	for(int c=0; c<neighborGrid.gridCells.size(); c++)
	{
		if(tmpIndex == -1)
		{
			minRow = neighborGrid.gridCells[c].rowNumber;
			minCol = neighborGrid.gridCells[c].colNumber;
			tmpIndex = c;
		}
		else if(neighborGrid.gridCells[c].rowNumber <= minRow && neighborGrid.gridCells[c].colNumber <= minCol)
		{
			minRow = neighborGrid.gridCells[c].rowNumber;
			minCol = neighborGrid.gridCells[c].colNumber;
			tmpIndex = c;
		}
	}
	
	for(int i=0; i<inliers.size(); i++)
	{
		Vec2f featurePos = aPts[inliers[i]];
		Vec2f a = featurePos*scale + translation;
		
		nTmpIndex = -1;
		
		//find closest grid center
		/*for(int c=0; c<neighborGrid.gridCells.size(); c++)
		{
			if(tmpIndex == -1)
			{
				tmpIndex = c;
				tmpDist = (neighborGrid.gridCells[c].center-featurePos).length();
			}
			else if((neighborGrid.gridCells[c].center-featurePos).length() < tmpDist)
			{
				tmpIndex = c;
				tmpDist = (neighborGrid.gridCells[c].center-featurePos).length();
			}
		}*/
		
		Vec2f diff = neighborGrid.gridCells[tmpIndex].center-featurePos;
		
		diff = diff*scale;
		
		Vec2f b = bPts[inliers[i]];
		featurePos = b + diff;
		
		//find closest grid center
		for(int c=0; c<refGrid.gridCells.size(); c++)
		{
			if(nTmpIndex == -1)
			{
				nTmpIndex = c;
				nTmpDist = (refGrid.gridCells[c].center-featurePos).length();
			}
			else if((refGrid.gridCells[c].center-featurePos).length() < nTmpDist)
			{
				nTmpIndex = c;
				nTmpDist = (refGrid.gridCells[c].center-featurePos).length();
			}
		}
		
		//printf("dist:%f, scale:%f\n", nTmpDist, scale);
		
		if(nTmpDist < repThresh)
		{
			rowShift = refGrid.gridCells[nTmpIndex].rowNumber - neighborGrid.gridCells[tmpIndex].rowNumber;
			colShift = refGrid.gridCells[nTmpIndex].colNumber - neighborGrid.gridCells[tmpIndex].colNumber;
			
			//if(rowShift > (refGrid.numRows - 1) || rowShift < (1 - neighborGrid.numRows))
			//	continue;
			//if(colShift > (refGrid.numColumns - 1) || colShift < (1 - neighborGrid.numColumns))
			//	continue;
			
			int iter;
			for(iter=0; iter<alignments.size(); iter++)
			{
				if(alignments[iter][0] == rowShift && alignments[iter][1] == colShift)
				{
					alignmentVotes[iter] += 1;
				}
			}
			if(iter == alignments.size())
			{
				alignments.push_back(Vec2i(rowShift, colShift));
				alignmentVotes.push_back(1);
			}
		}
		else if(nTmpDist < width && nTmpDist < height)
		{
			//the corresponding cell could be missing
			Vec2f pos = refGrid.gridCells[nTmpIndex].center;
			int noRowsUp = 0;
			int noRowsDown = 0;
			if(refGrid.yTrans[1] != 0.0)
			{
				noRowsUp = floor((float)nTmpDist / (float)refGrid.yTrans[1] + 0.5);
				noRowsDown = floor((float)nTmpDist / (float)refGrid.yTrans[1] + 0.5);
			}
			int noColsLeft = 0;
			int noColsRight = 0;
			if(refGrid.xTrans[0] != 0.0)
			{
				noColsLeft = floor((float)nTmpDist / (float)refGrid.xTrans[0] + 0.5);
				noColsRight = floor((float)nTmpDist / (float)refGrid.xTrans[0] + 0.5);
			}
			int newMinRowIndex = refGrid.gridCells[nTmpIndex].rowNumber;
			int newMinColIndex = refGrid.gridCells[nTmpIndex].colNumber;
			float newMinDistance = nTmpDist;
			for(int y=-noRowsUp; y<=noRowsDown; y++)
			{
				for(int x=-noColsLeft; x<=noColsRight; x++)
				{
					Vec2f newPos = pos + refGrid.yTrans*y + refGrid.xTrans*x;
					if((newPos-featurePos).length() < newMinDistance)
					{
						newMinDistance = (newPos-featurePos).length();
						newMinRowIndex = refGrid.gridCells[nTmpIndex].rowNumber + y;
						newMinColIndex = refGrid.gridCells[nTmpIndex].colNumber + x;
					}
				}
			}
			if(newMinDistance < repThresh)
			{
				rowShift = newMinRowIndex - neighborGrid.gridCells[tmpIndex].rowNumber;
				colShift = newMinColIndex - neighborGrid.gridCells[tmpIndex].colNumber;
				
				//if(rowShift > (refGrid.numRows - 1) || rowShift < (1 - neighborGrid.numRows))
				//	continue;
				//if(colShift > (refGrid.numColumns - 1) || colShift < (1 - neighborGrid.numColumns))
				//	continue;
				int iter;
				for(iter=0; iter<alignments.size(); iter++)
				{
					if(alignments[iter][0] == rowShift && alignments[iter][1] == colShift)
					{
						alignmentVotes[iter] += 1;
					}
				}
				if(iter == alignments.size())
				{
					alignments.push_back(Vec2i(rowShift, colShift));
					alignmentVotes.push_back(1);
				}
			}

		}
	}
	
	int maxVotes = 0;
	int maxVoteIndex = -1;
	for(int i=0; i<alignments.size(); i++)
	{
		//printf("alignment (%d %d) has vote %d\n", alignments[i][0], alignments[i][1], alignmentVotes[i]);
		if(maxVoteIndex == -1)
		{
			maxVoteIndex = i;
			maxVotes = alignmentVotes[i];
			rowShift = alignments[i][0];
			colShift = alignments[i][1];
		}
		else if(alignmentVotes[i] > maxVotes)
		{
			maxVoteIndex = i;
			maxVotes = alignmentVotes[i];
			rowShift = alignments[i][0];
			colShift = alignments[i][1];
		}
	}
	
	if(maxVoteIndex != -1)
		return maxVotes;
	else
		return -1;
}

float RepetitionFinder::translateTransformationToGridRANSAC(repetitionGrid &refGrid, repetitionGrid &neighborGrid, vector<Vec2f> &aPts, vector<Vec2f> &bPts, vector<float> &scale, vector<Vec2f> &translation, 
															vector<vector<int> > &inlierMatches, float repThresh, int width, int height, vector<int> &shiftCol, vector<int> &shiftRow, vector<vector<int> > &bestMatches, 
															int &bestCandidate)
{
	int totalMatches = aPts.size();
	
	vector<Vec2i> shifts;
	vector<int> scores;
	vector<int> scoreIndices;
	vector<int> totalInliers;
	
	int rowShift, colShift;
	
	//find shift defined for scale/translation
	for(int i=0; i<scale.size(); i++)
	{
		int success = translateTransformationToGrid(refGrid, neighborGrid, aPts, bPts, inlierMatches[i], scale[i], translation[i], repThresh, width, height, colShift, rowShift);
		if(success != -1)
		{
			int iter;
			for(iter=0; iter<shifts.size(); iter++)
			{
				if(shifts[iter][0]==colShift && shifts[iter][1]==rowShift)
				{
					if(success > scores[iter])
					{
						scores[iter] = success;
						scoreIndices[iter] = i;
					}
					break;
				}
			}
			if(iter == shifts.size())
			{
				shifts.push_back(Vec2i(colShift, rowShift));
				scores.push_back(success);
				scoreIndices.push_back(i);
			}
		}
	}
	
	int maxCount = 0;
	int maxCountIndex = -1;
	
	for(int i=0; i<scores.size(); i++)
	{
		printf("shift row:%d, shift col:%d, count:%d, ratio:%f\n", shifts[i][1], shifts[i][0], scores[i], (float)scores[i] / (float)aPts.size());
		int inlierCount;
		int ind;
		if(scores[i] > 12)
		{
			vector<int>::iterator iterRow = find(shiftRow.begin(), shiftRow.end(), shifts[i][1]);
			vector<int>::iterator iterCol = find(shiftCol.begin(), shiftCol.end(), shifts[i][0]);
			
			if(iterRow != shiftRow.end() && iterCol != shiftCol.end() && (iterRow-shiftRow.begin()) == (iterCol - shiftCol.begin()))
			{
				ind = iterRow - shiftRow.begin();
				for(int j=0; j<inlierMatches[scoreIndices[i]].size(); j++)
				{
					if(find(bestMatches[ind].begin(), bestMatches[ind].end(), inlierMatches[scoreIndices[i]][j]) == bestMatches[ind].end())
					{
						bestMatches[ind].push_back(inlierMatches[scoreIndices[i]][j]);
					}
				}
				inlierCount = bestMatches[ind].size();
			}
			else
			{
				shiftRow.push_back(shifts[i][1]);
				shiftCol.push_back(shifts[i][0]);
				
				vector<int> bm;
				for(int j=0; j<inlierMatches[scoreIndices[i]].size(); j++)
				{
					if(find(bm.begin(), bm.end(), inlierMatches[scoreIndices[i]][j]) == bm.end())
					{
						bm.push_back(inlierMatches[scoreIndices[i]][j]);
					}
				}
				bestMatches.push_back(bm);
				inlierCount = bm.size();
				ind = bestMatches.size() - 1;
			}
			if(maxCountIndex == -1)
			{
				maxCountIndex = ind;
				maxCount = inlierCount;
			}
			else if(inlierCount > maxCount)
			{
				maxCountIndex = ind;
				maxCount = inlierCount;
			}
		}
	}
	
	bestCandidate = maxCountIndex;
	
	float ratio = (float)(maxCount) / (float)(totalMatches);
	
	if(maxCount < 20)
		ratio = -1.0;
	
	return ratio;
}

bool RepetitionFinder::refineMatching(repetitionGrid &refGrid, repetitionGrid &neighborGrid, repetitionGrid &currentGrid, Vec2i neighborShift, float repThresh, 
									  float scale, Vec2f &translate, float ratio, vector<int> &shiftCol, vector<int> &shiftRow)
{
	int numRefRows = currentGrid.numRows;
	int numRefCols = currentGrid.numColumns;
	vector<vector<int> > refDescriptors;
	refDescriptors.resize(numRefRows);
	for(int i=0; i<numRefRows; i++)
	{
		refDescriptors[i].resize(numRefCols, -1);
	}
	for(int i=0; i<currentGrid.gridCells.size(); i++)
	{
		refDescriptors[currentGrid.gridCells[i].rowNumber][currentGrid.gridCells[i].colNumber] = i;//currentGrid.gridCells[i].descriptor;
	}
	
	int numNRows = refGrid.numRows;
	int numNCols = refGrid.numColumns;
	vector<vector<float> > nDescriptors;
	nDescriptors.resize(numNRows);
	for(int i=0; i<numNRows; i++)
	{
		nDescriptors[i].resize(numNCols, -1);
	}
	for(int i=0; i<refGrid.gridCells.size(); i++)
	{
		nDescriptors[refGrid.gridCells[i].rowNumber][refGrid.gridCells[i].colNumber] = i;//refGrid.gridCells[i].descriptor;
	}

	//find shift defined for scale/translation
	int rowShift = 0; int colShift = 0;
	int refIndex = -1; int closestIndex = -1;
	float dist = repThresh;
	for(int c=0; c<refGrid.gridCells.size(); c++)
	{
		Vec2f pos = refGrid.gridCells[c].center*scale + translate;
		
		for(int n=0; n<neighborGrid.gridCells.size(); n++)
		{
			if(closestIndex == -1)
			{
				closestIndex = n;
				refIndex = c;
				dist = (pos-neighborGrid.gridCells[n].center).length();
			}
			else if((pos-neighborGrid.gridCells[n].center).length() < dist)
			{
				closestIndex = n;
				refIndex = c;
				dist = (pos-neighborGrid.gridCells[n].center).length();
			}
		}
	}
	
	int minCol, maxCol, minRow, maxRow;
	//search in local neighborhood
	int colRange = (numRefCols + numNCols - 1) / 3;
	int rowRange = (numRefRows + numNRows - 1) / 3;
	
	if(dist < repThresh)
	{
		rowShift = neighborGrid.gridCells[closestIndex].rowNumber - refGrid.gridCells[refIndex].rowNumber;
		colShift = neighborGrid.gridCells[closestIndex].colNumber - refGrid.gridCells[refIndex].colNumber;
		rowShift += neighborShift[0];
		colShift += neighborShift[1];
		printf("initial rowShift:%d colShift:%d\n", rowShift, colShift);
		minCol = max(colShift-colRange, -(numNCols-1));
		maxCol = min(colShift+colRange, numRefCols);
		minRow = max(rowShift-rowRange, -(numNRows-1));
		maxRow = min(rowShift+rowRange, numRefRows);
	}
	else
	{
		printf("Not enough repetitions to refine the matching.\n");
		minCol = -(numNCols-1); maxCol = numRefCols;
		minRow = -(numNRows-1); maxRow = numRefRows;
	}
	
	if(ratio > 0.6)
	{
		shiftCol.push_back(colShift);
		shiftRow.push_back(rowShift);
		return true;
	}
	else if(ratio < 0.4)
	{
		minCol = -(numNCols-1); maxCol = numRefCols;
		minRow = -(numNRows-1); maxRow = numRefRows;
	}
	
	vector<vector<float> > scores;
	scores.resize(numRefRows+(numNRows-1));
	for(int i=0; i<scores.size(); i++)
		scores[i].resize(numRefCols+(numNCols-1), FLT_MAX);
	
	int countThreshold = floor(min(numRefCols, numNCols)/2.0)*floor(min(numRefRows, numNRows)/2.0);
	
	for(int x=minCol; x<maxCol; x++)
	{
		for(int y=minRow; y<maxRow; y++)
		{
			//shift x cols and y rows
			scores[y+(numNRows-1)][x+(numNCols-1)] = 0;
			//printf("shifting %d cols and %d rows\n", x, y);
			int count = 0;
			for(int i=0; i<numNCols; i++)
			{
				for(int j=0; j<numNRows; j++)
				{
					if(x+i >=0 && x+i < numRefCols && y+j >= 0 && y+j < numRefRows)
					{
						if(refDescriptors[y+j][x+i] != -1 && nDescriptors[j][i] != -1)
						{
							float diff = 0.0;
							for(int c=0; c<currentGrid.gridCells[refDescriptors[y+j][x+i]].pcaCoord.size(); c++)
							{
								diff += (currentGrid.gridCells[refDescriptors[y+j][x+i]].pcaCoord[c] - refGrid.gridCells[nDescriptors[j][i]].pcaCoord[c]) *
								(currentGrid.gridCells[refDescriptors[y+j][x+i]].pcaCoord[c] - refGrid.gridCells[nDescriptors[j][i]].pcaCoord[c]);
							}
							
							//printf("ref:%f n:%f\n", refDescriptors[y+j][x+i], nDescriptors[j][i]);
							scores[y+(numNRows-1)][x+(numNCols-1)] += diff;//(refDescriptors[y+j][x+i] - nDescriptors[j][i])*(refDescriptors[y+j][x+i] - nDescriptors[j][i]);
							if(refDescriptors[y+j][x+i] != -1 && nDescriptors[j][i] != -1)
								count ++;
						}
					}
					//else
					// {
					// scores[y+(numNRows-1)][x+(numNCols-1)] += nDescriptors[j][i] * nDescriptors[j][i];
					// }
				}
			}
			if(count < countThreshold)
				scores[y+(numNRows-1)][x+(numNCols-1)] = FLT_MAX;
			else
				scores[y+(numNRows-1)][x+(numNCols-1)] = scores[y+(numNRows-1)][x+(numNCols-1)]*float(refGrid.gridCells.size()*neighborGrid.gridCells.size())/float(count*count);
			//printf("score at col:%d, row:%d:%f with count:%d\n", x, y, scores[y+(numNRows-1)][x+(numNCols-1)], count);
		}
	}
	
	vector<gridRegistrationCandidate> candidates;
	
	for(int i=0; i<scores.size(); i++)
	{
		for(int j=0; j<scores[i].size(); j++)
		{
			if(scores[i][j] == countThreshold)
				continue;
			else
			{
				gridRegistrationCandidate can;
				can.colShift = j;
				can.rowShift = i;
				can.score = scores[i][j];
				candidates.push_back(can);
			}
		}
	}
	
	sort(candidates.begin(), candidates.end());
	
	int k=min(5, (int)(candidates.size()));
	
	for(int i=0; i<k; i++)
	{
		shiftCol.push_back(candidates[i].colShift - (numNCols-1));
		shiftRow.push_back(candidates[i].rowShift - (numNRows-1));
	}
	printf("max score at row %d and col %d:%f\n", shiftRow[0], shiftCol[0], candidates[0].score);
	return true;
}

void RepetitionFinder::findBestMatchingScore(repetitionGrid &refGrid, repetitionGrid &neighborGrid, float &score, float &ratio, int &shiftCol, int &shiftRow)
{
	int numRefRows = refGrid.numRows;
	int numRefCols = refGrid.numColumns;
	vector<vector<int> > refDescriptors;
	refDescriptors.resize(numRefRows);
	for(int i=0; i<numRefRows; i++)
	{
		refDescriptors[i].resize(numRefCols, -1);
	}
	for(int i=0; i<refGrid.gridCells.size(); i++)
	{
		refDescriptors[refGrid.gridCells[i].rowNumber][refGrid.gridCells[i].colNumber] = i;//refGrid.gridCells[i].descriptor;
	}
	
	int numNRows = neighborGrid.numRows;
	int numNCols = neighborGrid.numColumns;
	vector<vector<int> > nDescriptors;
	nDescriptors.resize(numNRows);
	for(int i=0; i<numNRows; i++)
	{
		nDescriptors[i].resize(numNCols, -1);
	}
	for(int i=0; i<neighborGrid.gridCells.size(); i++)
	{
		nDescriptors[neighborGrid.gridCells[i].rowNumber][neighborGrid.gridCells[i].colNumber] = i;//neighborGrid.gridCells[i].descriptor;
	}
	
	vector<vector<float> > scores;
	scores.resize(numRefRows+(numNRows-1));
	for(int i=0; i<scores.size(); i++)
		scores[i].resize(numRefCols+(numNCols-1), 0.0);
	
	int countThreshold = floor(min(numRefCols, numNCols)/2.0)*floor(min(numRefRows, numNRows)/2.0);
	
	for(int x=-(numNCols-1); x<numRefCols; x++)
	{
		for(int y=-(numNRows-1); y<numRefRows; y++)
		{
			//shift x cols and y rows
			//printf("shifting %d cols and %d rows\n", x, y);
			int count = 0;
			for(int i=0; i<numNCols; i++)
			{
				for(int j=0; j<numNRows; j++)
				{
					if(x+i >=0 && x+i < numRefCols && y+j >= 0 && y+j < numRefRows)
					{
						if(refDescriptors[y+j][x+i] != -1 && nDescriptors[j][i] != -1)
						{
							//printf("ref:%d n:%d\n", refDescriptors[y+j][x+i], nDescriptors[j][i]);
							float diff = 0.0;
							for(int c=0; c<refGrid.gridCells[refDescriptors[y+j][x+i]].pcaCoord.size(); c++)
							{
								diff += (refGrid.gridCells[refDescriptors[y+j][x+i]].pcaCoord[c] - neighborGrid.gridCells[nDescriptors[j][i]].pcaCoord[c]) *
										(refGrid.gridCells[refDescriptors[y+j][x+i]].pcaCoord[c] - neighborGrid.gridCells[nDescriptors[j][i]].pcaCoord[c]);
							}
							scores[y+(numNRows-1)][x+(numNCols-1)] += diff;//(refDescriptors[y+j][x+i] - nDescriptors[j][i])*(refDescriptors[y+j][x+i] - nDescriptors[j][i]);
							if(refDescriptors[y+j][x+i] != -1 && nDescriptors[j][i] != -1)
								count ++;
						}
						else 
						{
							//printf("refDesc:%d, nDesc:%d\n", refDescriptors[y+j][x+i], nDescriptors[j][i]);
							int deneme = 1;
						}

					}
					else
					{
						//printf("n:%f\n", nDescriptors[j][i]);
						//scores[y+(numNRows-1)][x+(numNCols-1)] += nDescriptors[j][i] * nDescriptors[j][i];
						//printf("ref:%d n:%d\n", refDescriptors[y+j][x+i], nDescriptors[j][i]);
						float diff = 0.0;
						for(int c=0; c<neighborGrid.gridCells[nDescriptors[j][i]].pcaCoord.size(); c++)
						{
							diff += (neighborGrid.gridCells[nDescriptors[j][i]].pcaCoord[c]) *
									(neighborGrid.gridCells[nDescriptors[j][i]].pcaCoord[c]);
						}
						scores[y+(numNRows-1)][x+(numNCols-1)] += 75.0*75.0;//diff;
					}
				}
			}
			if(count < countThreshold)
				scores[y+(numNRows-1)][x+(numNCols-1)] = FLT_MAX;//countThreshold;//numRefCols*numRefRows; //arbitrary high score
			else
				scores[y+(numNRows-1)][x+(numNCols-1)] = scores[y+(numNRows-1)][x+(numNCols-1)];//*float(refGrid.gridCells.size()*neighborGrid.gridCells.size())/float(count*count);
			//printf("score at col:%d, row:%d:%f with count:%d\n", x, y, scores[y+(numNRows-1)][x+(numNCols-1)], count);
		}
	}
	
	float maxScore, secondMaxScore;
	int maxScoreRow, maxScoreCol;
	
	maxScore = -1.0; secondMaxScore = -1.0;
	
	for(int i=0; i<scores.size(); i++)
	{
		for(int j=0; j<scores[i].size(); j++)
		{
			if(scores[i][j] == countThreshold)
				continue;
			
			if(maxScore == -1.0)
			{
				maxScore = scores[i][j];
				maxScoreRow = i;
				maxScoreCol = j;
			}
			else if(secondMaxScore == -1.0)
			{
				if(scores[i][j] < maxScore)
				{
					secondMaxScore = maxScore;
					maxScore = scores[i][j];
					maxScoreRow = i;
					maxScoreCol = j;
				}
				else
				{
					secondMaxScore = scores[i][j];
				}
			}
			else
			{
				if(scores[i][j] < maxScore)
				{
					secondMaxScore = maxScore;
					maxScore = scores[i][j];
					maxScoreRow = i;
					maxScoreCol = j;
				}
				else if(scores[i][j] < secondMaxScore)
				{
					secondMaxScore = scores[i][j];
				}
			}
		}
	}
	
	printf("max score at row %d and col %d:%f, second max score:%f, ratio:%f\n", maxScoreRow-(numNRows-1), maxScoreCol-(numNCols-1), maxScore, secondMaxScore, secondMaxScore/maxScore);
	shiftCol = maxScoreCol - (numNCols-1);
	shiftRow = maxScoreRow - (numNRows-1);
	score = maxScore;
	ratio = secondMaxScore/maxScore;
}

repetitionGrid RepetitionFinder::mergeGrids(repetitionGrid &grid1, repetitionGrid &grid2, int shiftCol, int shiftRow)
{
	repetitionGrid currentGrid;
	
	if(shiftCol < 0)
	{
		currentGrid.numColumns = grid1.numColumns + abs(shiftCol);
		if(grid2.numColumns - abs(shiftCol) > grid1.numColumns)
		{
			currentGrid.numColumns += grid2.numColumns - abs(shiftCol) - grid1.numColumns;
		}
	}
	else
	{
		currentGrid.numColumns = grid1.numColumns;
		if(grid2.numColumns + shiftCol > grid1.numColumns)
		{
			currentGrid.numColumns += (grid2.numColumns + shiftCol - grid1.numColumns);
		}
	}
	
	if(shiftRow < 0)
	{
		currentGrid.numRows = grid1.numRows + abs(shiftRow);
		if(grid2.numRows - abs(shiftRow) > grid1.numRows)
		{
			currentGrid.numRows += grid2.numRows - abs(shiftRow) - grid1.numRows;
		}
	}
	else
	{
		currentGrid.numRows = grid1.numRows;
		if(grid2.numRows + shiftRow > grid1.numRows)
		{
			currentGrid.numRows += (grid2.numRows + shiftRow - grid1.numRows);
		}
	}
	
	vector<vector<int> > visImages;
	visImages.resize(currentGrid.numRows);
	
	for(int i=0; i<visImages.size(); i++)
	{
		visImages[i].resize(currentGrid.numColumns, 0);
	}
	
	int c = min(0, shiftCol);
	int r = min(0, shiftRow);
	
	c = abs(c); r = abs(r);
	
	for(int i=0; i<grid1.gridCells.size(); i++)
	{
		visImages[grid1.gridCells[i].rowNumber+r][grid1.gridCells[i].colNumber+c] += grid1.gridCells[i].visImages;
	}
	
	c = max(0, shiftCol);
	r = max(0, shiftRow);
	
	for(int i=0; i<grid2.gridCells.size(); i++)
	{
		visImages[grid2.gridCells[i].rowNumber+r][grid2.gridCells[i].colNumber+c] += grid2.gridCells[i].visImages;
	}
	
	float countTotal = 0;
	for(int i=0; i<visImages.size(); i++)
	{
		for(int j=0; j<visImages[i].size(); j++)
		{
			if(visImages[i][j] != 0)
			{
				currentGrid.gridCells.push_back(repetitionCell(Vec2f(0.0, 0.0)));
				currentGrid.gridCells[currentGrid.gridCells.size()-1].rowNumber = i;
				currentGrid.gridCells[currentGrid.gridCells.size()-1].colNumber = j;
				currentGrid.gridCells[currentGrid.gridCells.size()-1].visImages = visImages[i][j];
				countTotal += 1.0;
			}
		}
	}
	
	/*Img tmp(currentGrid.numColumns*100, currentGrid.numRows*100);
	for(int i=0; i<currentGrid.numRows; i++)
	{
		for(int j=0; j<currentGrid.numColumns; j++)
		{
			Vec3uc c;
			if(scores[i][j] == 0.0)
				c = Color(0.0, 0.0, 0.0);
			else
				c = MathUtils::generateColorFromValue(scores[i][j], 0.0, 1.0);
			for(int y=i*100; y<(i+1)*100; y++)
			{
				for(int x=j*100; x<(j+1)*100; x++)
				{
					tmp.setColor(x, y, Color(c[0]/255.0, c[1]/255.0, c[2]/255.0));
				}
			}
		}
	}
	tmp.write("/Users/ceylan/Desktop/current.png");*/
	return currentGrid;
}
