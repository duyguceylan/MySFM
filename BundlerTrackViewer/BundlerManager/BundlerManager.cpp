/*
 *  BundlerManager.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/15/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include <sys/stat.h>
#include "BundlerManager.h"
#include "BundlerApp.h"
#include "Image.h"
#include "../Parser/BundlerOutputParser.h"
#include "../MathUtils/GeometryUtils.h"
#include "../MathUtils/MathUtils.h"
#include "../TextureUtils/TextureManager.h"
#include "SIFT.h"
#include "MRFDepthLabeling.h"
#include "opencv2/opencv.hpp"
#include <exiv2/exiv2.hpp>

#define KODAK_CCD_WIDTH 6.16
#define IPHONE4S_CCD_WIDTH 4.54

#define rectImgWidth 2000;
#define rectImgHeight 1500;

BundlerManager::BundlerManager(string _mainDirectory)
{
	mainDirectory = _mainDirectory;
	bundlerOutputPath = mainDirectory + "/comparisonJiang2012/bundle.out";
	
	keySpaceAllocated = false;
	lineDetector = NULL;
	kpExtractor = NULL;
	kpMatcher = NULL;
	geoManager = NULL;
	boParser = NULL;
	imgRect = NULL;
	
	lines.clear();
	keys.clear();
	keyInfo.clear();
	noFeaturePoints.clear();
	
	points3D.clear();
	
	gridStrokes.clear();
	
	string denseFilePath = mainDirectory + "/models/dense.ply";
	ifstream finDense(denseFilePath.c_str(), ios::in);
	if(finDense)
	{
		densePC = new PointCloud(denseFilePath.c_str());
	}
	
	string sparseFilePath = mainDirectory + "/models/sparse.ply";
	ifstream finSparse(sparseFilePath.c_str(), ios::in);
	if(finSparse)
	{
		sparsePC = new PointCloud(sparseFilePath.c_str());
	}
	
	printf("main directory:%s\n", mainDirectory.c_str());
	
	boParser = new BundlerOutputParser();
	//boParser->convertSfMOutputToPMVS("/Users/ceylan/MVSReconstruction/SfM/SfM/data/museum/model-0-cams.txt");
	
}

void BundlerManager::addNewImagePath(std::string path)
{
	imageFileNames.push_back(path); 
	int width, height;
	float focalLength;
	
	Exiv2::Image::AutoPtr tmp = Exiv2::ImageFactory::open(path.c_str());
    tmp->readMetadata();
	
    Exiv2::ExifData &exifData = tmp->exifData();
    if (exifData.empty()) 
	{
        std::string error(path.c_str());
        error += ": No Exif data found in the file";
		return ;
    }
    Exiv2::ExifData::const_iterator end = exifData.end();
    for (Exiv2::ExifData::const_iterator i = exifData.begin(); i != end; ++i) 
	{
		string attr(i->key());
		if(attr.compare("Exif.Photo.PixelXDimension") == 0)
		{
			width = (int)(i->value().toFloat());
		}
		else if(attr.compare("Exif.Photo.PixelYDimension") == 0)
		{
			height = (int)(i->value().toFloat());
		}
		else if(attr.compare("Exif.Photo.FocalLength") == 0)
		{
			focalLength = i->value().toFloat();
		}
    }
	 
	
	printf("width:%d, height:%d, focalLength:%f\n", width, height, focalLength);
	
	int maxSize = width;
	if(height > width)
		maxSize = height;
	
	Matrix3f K;
	K.setIdentity();
	float focalLengthInPixels = focalLength*(float)maxSize/KODAK_CCD_WIDTH;
	K[0][0] = focalLengthInPixels;
	K[1][1] = focalLengthInPixels;
	K[2][0] = (float)(width)/2.0;
	K[2][1] = (float)(height)/2.0;
	
	intrinsicMatrices.push_back(K);
	
	imageSizes.push_back(Vec2i(width, height));
	gridStrokes.push_back(vector<vector<Vec2f> >());
}

void BundlerManager::computeEdges()
{
	if(lineDetector == NULL)
	{
		lineDetector = new LineDetector();
	}
	for(int i=0; i<lines.size(); i++)
	{
		lines[i].clear();
	}
	lines.clear();
	lines.resize(imageFileNames.size());
	
	for(int i=0; i<imageFileNames.size(); i++)
	{
		Photo p(imageFileNames[i]);
		float sigma = sqrt(2.0);
		float tlow = 0.07;
		float thigh = 0.1;
		float lengthThreshold = 30;
		
		char lineFile[1024];
		sprintf(lineFile, "%08d.%s", i, "txt");
		string filename = mainDirectory + "/images/" + lineFile;
		
		ifstream fin(filename.c_str(), ios::in);
		if(fin)
		{
			lineDetector->read2DLineSegments(filename.c_str(), lines[i], lengthThreshold);
		}
		else
		{
			lineDetector->find2DLineSegments(&p, sigma, tlow, thigh, lengthThreshold, lines[i], filename.c_str());
		}
	}
}

Img* BundlerManager::RectifyImage(int index, int& numberOfPlanes)
{
	if(imgRect == NULL)
	{
		imgRect = new ImageRectification();
	}
	if(rectifiedImages.size() == 0)
	{
		rectifiedImages.resize(imageFileNames.size());
		noRectifiedImages.resize(imageFileNames.size(), 0);
		rectifyingHomographies.resize(imageFileNames.size());
		vanishingPointDist.resize(imageFileNames.size(),Vec2f(0.0,0.0));
		
		for(int i=0; i<imageFileNames.size(); i++)
		{
			rectifiedImages[i].resize(2, NULL);
			rectifyingHomographies[i].resize(2);
		}
	}
	
	if(index > imageFileNames.size() -1)
	{
		printf("Image index out of bounds.\n");
		return NULL;
	}
	
	if(rectifiedImages[index][0] == NULL)
	{
		char imgFile[1024];
		sprintf(imgFile, "%08d_0.%s", index, "png");
		string filename = mainDirectory + "/rectification/" + imgFile;
		char homFile[1024];
		sprintf(homFile, "homography%d_0.txt", index);
		string homographyFileName = mainDirectory + "/rectification/" + homFile;
		
		char newImgFile[1024];
		sprintf(newImgFile, "%08d_1.%s", index, "png");
		string newFilename = mainDirectory + "/rectification/" + newImgFile;
		char newHomFile[1024];
		sprintf(newHomFile, "homography%d_1.txt", index);
		string newHomographyFileName = mainDirectory + "/rectification/" + newHomFile;
		
		ifstream fin(filename.c_str());
		if(fin)
		{
			rectifiedImages[index][0] = new Img();
			rectifiedImages[index][0]->read(filename.c_str());
			Matrix3f H;
			ifstream homFile(homographyFileName.c_str(), ios::in);
			homFile >> H[0][0] >> H[1][0] >> H[2][0];
			homFile >> H[0][1] >> H[1][1] >> H[2][1];
			homFile >> H[0][2] >> H[1][2] >> H[2][2];
			homFile.close();
			rectifyingHomographies[index][0] = H;
			noRectifiedImages[index] = noRectifiedImages[index]+1;
			
			ifstream finSecond(newFilename.c_str());
			if(finSecond)
			{
				rectifiedImages[index][1] = new Img();
				rectifiedImages[index][1]->read(newFilename.c_str());
				Matrix3f HSecond;
				ifstream newHomFile(newHomographyFileName.c_str(), ios::in);
				newHomFile >> HSecond[0][0] >> HSecond[1][0] >> HSecond[2][0];
				newHomFile >> HSecond[0][1] >> HSecond[1][1] >> HSecond[2][1];
				newHomFile >> HSecond[0][2] >> HSecond[1][2] >> HSecond[2][2];
				newHomFile.close();
				rectifyingHomographies[index][1] = HSecond;
				noRectifiedImages[index] = noRectifiedImages[index]+1;
			}
		}
		else
		{
			Img im;
			im.read(imageFileNames[index].c_str());
			imgRect->setImage(&im);
			bool secondHorizontalVPFound = false;
			Matrix3f H = imgRect->rectifyImage(&(rectifiedImages[index][0]), vanishingPointDist[index], secondHorizontalVPFound);
			rectifyingHomographies[index][0] = H;
		
			//save rectified image
			rectifiedImages[index][0]->write(filename.c_str());
			//save homography
			ofstream fout(homographyFileName.c_str(), ios::out);
			fout << H[0][0] << " " << H[1][0] << " " << H[2][0] << endl;
			fout << H[0][1] << " " << H[1][1] << " " << H[2][1] << endl;
			fout << H[0][2] << " " << H[1][2] << " " << H[2][2] << endl;
			fout.close();
			noRectifiedImages[index] = noRectifiedImages[index]+1;
			
			if(secondHorizontalVPFound)
			{
				Vec2f vanishingPtDistance;
				Matrix3f newRectMatrix = imgRect->rectifyImageWrtSecondHorizontalDirection(&(rectifiedImages[index][1]), vanishingPtDistance);
				
				rectifiedImages[index][1]->write(newFilename.c_str());
				rectifyingHomographies[index][1] = newRectMatrix;
				
				ofstream foutHom(newHomographyFileName.c_str(), ios::out);
				foutHom << newRectMatrix[0][0] << " " << newRectMatrix[1][0] << " " << newRectMatrix[2][0] << endl;
				foutHom << newRectMatrix[0][1] << " " << newRectMatrix[1][1] << " " << newRectMatrix[2][1] << endl;
				foutHom << newRectMatrix[0][2] << " " << newRectMatrix[1][2] << " " << newRectMatrix[2][2] << endl;
				foutHom.close();
				
				noRectifiedImages[index] = noRectifiedImages[index]+1;
			}
		}
	}
	numberOfPlanes = noRectifiedImages[index];
	return rectifiedImages[index][0];
}

Img* BundlerManager::getRectifiedImage(int imgIndex, int planeIndex)
{
	if(planeIndex >= noRectifiedImages[imgIndex])
		return NULL;
	return rectifiedImages[imgIndex][planeIndex];
}

void BundlerManager::readGridMatches(int index)
{
	char imgFile[1024];
	sprintf(imgFile, "%08d.%s", index, "txt");
	string filename = mainDirectory + "/matches/" + imgFile;
	
	RepetitionFinder rf;
	
	
	ifstream fin(filename.c_str(), ios::in);
	if(fin)
	{
		int noTemplates;
		int planeId;
		int noGrids;
		int noMatches;
		int w, h;
		float x, y, score;
		fin >> noTemplates;
		if(repetitions.size() == 0)
		{
			repetitions.resize(noTemplates);
			refinedRepetitions.resize(noTemplates);
			repeatingGrids.resize(noTemplates);
			elementSizes.resize(noTemplates);
			
			for(int t=0; t<noTemplates; t++)
			{
				int noImages = imageFileNames.size();
				elementSizes[t].resize(noImages);
				repetitions[t].resize(noImages);
				refinedRepetitions[t].resize(noImages);
				repeatingGrids[t].resize(noImages);
			}
		}
		
		for(int t=0; t<noTemplates; t++)
		{
			elementSizes[t][index].resize(noRectifiedImages[index]);
			repetitions[t][index].resize(noRectifiedImages[index]);
			refinedRepetitions[t][index].resize(noRectifiedImages[index]);
			repeatingGrids[t][index].resize(noRectifiedImages[index]);
			
			for(int p=0; p<noRectifiedImages[index]; p++)
			{
				fin >> noMatches >> w >> h;
				elementSizes[t][index][p] = Vec2i(w,h);
				if(noMatches == 0)
				{
					continue;
				}
		
				vector<Vec2f> tempMatches;
				vector<float> scores;
			
				for(int i=0; i<noMatches; i++)
				{
					fin	>> x >> y >> score;
					tempMatches.push_back(Vec2f(x,y));
					scores.push_back(score);
				}
			
				vector<vector<Vec2f> > finalGroups;
				vector<Vec2f> finalTransformations;
			
				rf.findGridGeneratingVectors(tempMatches, scores, w/2.0, h/2.0, finalGroups, finalTransformations);
			
				repetitions[t][index][p].resize(finalGroups.size());
				refinedRepetitions[t][index][p].resize(finalGroups.size());
				
				for(int group=0; group<finalGroups.size(); group++)
				{
					repetitionGrid grid;
					grid.xTrans = Vec2f(finalTransformations[group][0], 0.0);
					grid.yTrans = Vec2f(0.0, finalTransformations[group][1]);
											
					vector<float> groupScores;
					groupScores.resize(finalGroups[group].size(), float(group+1)/(float)(finalGroups.size()));
					rf.segmentIntoRowsAndColumns(finalGroups[group], groupScores, grid, grid.yTrans, grid.xTrans, w, h);
					repeatingGrids[t][index][p].push_back(grid);
				
					for(int j=0; j<grid.gridCells.size(); j++)
					{
						Vec2f center = grid.gridCells[j].center;
					
						vector<Vec2f> coord;
						Vec3f pos;
						pos = rectifyingHomographies[index][p] * Vec3f(center[0] - w/2, center[1] - h/2, 1.0);
						coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
						pos = rectifyingHomographies[index][p] * Vec3f(center[0] + w/2, center[1] - h/2, 1.0);
						coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
						pos = rectifyingHomographies[index][p] * Vec3f(center[0] + w/2, center[1] + h/2, 1.0);
						coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
						pos = rectifyingHomographies[index][p] * Vec3f(center[0] - w/2, center[1] + h/2, 1.0);
						coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
					
						repetitions[t][index][p][group].push_back(coord);
				
						repeatingGrids[t][index][p][repeatingGrids[t][index][p].size()-1].gridCells[j].size = Vec2i(w,h);
					}
				}
			}
		}
	}			
}

void BundlerManager::readAllGridMatches()
{
	RepetitionFinder rf;
	
	int noImages = imageFileNames.size();
	
	repetitions.clear();
	repeatingGrids.clear();
	refinedRepetitions.clear();
	elementSizes.clear();
	
	for(int i=0; i<noImages; i++)
	{
		int noPlanes;
		RectifyImage(i, noPlanes);
		readGridMatches(i);
	}
}

void BundlerManager::findScaleOrder(int index, int planeIndex, vector<float> &scales)
{
	int noImages = imageFileNames.size();
	int noPlanes = 0;
	
	vector<vector<int> > planeIds;
	planeIds.resize(noImages);
	for(int i=0; i<noImages; i++)
	{
		for(int j=0; j<noRectifiedImages[i]; j++)
			planeIds[i].push_back(noPlanes+j);
		noPlanes += noRectifiedImages[i];
	}
	
	vector<int> processed;
	
	processed.push_back(planeIds[index][planeIndex]);
	scales.resize(noPlanes, -1);
	scales[planeIds[index][planeIndex]] = 1.0;
	bool stop = false;
	int maxCandidates = 29;
	while(!stop)
	{
		stop = true;
		
		vector<scaleEstimate> candidates;
		for(int p=0; p<processed.size(); p++)
		{
			
			
			int img, plane;
			bool found = false;
			for(int k=0; k<planeIds.size(); k++)
			{
				for(int m=0; m<planeIds[k].size(); m++)
				{
					if(processed[p] == planeIds[k][m])
					{
						img = k;
						plane = m;
						found = true;
						break;
					}
				}
				if(found)
					break;
			}
			
			bool changed = false;
			for(int c=0; c<maxCandidates; c++)
			{
				for(int i=0; i<noImages; i++)
				{
					for(int j=0; j<noRectifiedImages[i]; j++)
					{
						if(find(processed.begin(), processed.end(), planeIds[i][j]) != processed.end())
							continue;
						int curPlane = planeIds[i][j];
						if(scaleEstimates[curPlane].size() <= c)
							continue;
						if(scaleEstimates[curPlane][c].imgIndex2 == img && scaleEstimates[curPlane][c].planeIndex2 == plane)
						{
							candidates.push_back(scaleEstimates[curPlane][c]);
						}
					}
				}
			}
		}
		if(candidates.size() > 0)
		{
			stop = false;
			sort(candidates.begin(), candidates.end());
			processed.push_back(planeIds[candidates[0].imgIndex1][candidates[0].planeIndex1]);
			scales[planeIds[candidates[0].imgIndex1][candidates[0].planeIndex1]] = scales[planeIds[candidates[0].imgIndex2][candidates[0].planeIndex2]] * candidates[0].scale;
			printf("add img %d plane %d as a match of img %d plae %d with %d inliers\n", candidates[0].imgIndex1, candidates[0].planeIndex1, candidates[0].imgIndex2, candidates[0].planeIndex2, candidates[0].noInliers);
			printf("add to processed:%d with scale %f\n", processed[processed.size()-1], scales[planeIds[candidates[0].imgIndex1][candidates[0].planeIndex1]]);
		}
	}
}

void BundlerManager::estimateScales()
{
	SIFT *sift = new SIFT();
	RepetitionFinder rf;
	
	double threshold = 0.05;
	
	int noPlanes = 0;
	int noImages = imageFileNames.size();
	
	vector<vector<int> > planeIds;
	planeIds.resize(noImages);
	for(int i=0; i<noImages; i++)
	{
		for(int j=0; j<noRectifiedImages[i]; j++)
			planeIds[i].push_back(noPlanes+j);
		noPlanes += noRectifiedImages[i];
	}
	
	string scaleFilename = mainDirectory + "/matches/scales.txt";
	ifstream fin(scaleFilename.c_str());
	if(fin)
	{
		fin >> noPlanes;
		scaleEstimates.resize(noPlanes);
		int count = 0;
		for(int i=0; i<noImages; i++)
		{
			for(int p1=0; p1<noRectifiedImages[i]; p1++)
			{
				int noCandidates;
				int refImg, refPlane;
				fin >> refImg >> refPlane >> noCandidates;
				
				if(i != refImg || p1 != refPlane)
				{
					printf("Error in scales file!\n");
					return ;
				}
				
				for(int c=0; c<noCandidates; c++)
				{
					int img, plane, noInliers;
					float scale;
					fin >> img >> plane >> noInliers >> scale;
					scaleEstimate s;
					s.imgIndex1 = i; s.planeIndex1 = p1;
					s.imgIndex2 = img; s.planeIndex2 = plane;
					s.noInliers = noInliers; s.scale = scale;
					scaleEstimates[count].push_back(s);
				}
				count +=1 ;
			}
		}
		fin.close();
		
		for(int i=0; i<noImages; i++)
		{
			for(int p1=0; p1<noRectifiedImages[i]; p1++)
			{
				for(int j=i+1; j<noImages; j++)
				{
					for(int p2=0; p2<noRectifiedImages[j]; p2++)
					{
						int cInd, cInd2;
						
						for(int c=0; c<scaleEstimates[planeIds[i][p1]].size(); c++)
						{
							if(scaleEstimates[planeIds[i][p1]][c].imgIndex2 == j && scaleEstimates[planeIds[i][p1]][c].planeIndex2 == p2)
							{
								cInd = c;
							}
						}
						for(int c2=0; c2<scaleEstimates[planeIds[j][p2]].size(); c2++)
						{
							if(scaleEstimates[planeIds[j][p2]][c2].imgIndex2 == i && scaleEstimates[planeIds[j][p2]][c2].planeIndex2 == p1)
							{
								cInd2 = c2;
							}
						}
						
						if(cInd >= scaleEstimates[planeIds[i][p1]].size() || cInd2 >= scaleEstimates[planeIds[j][p2]].size())
							continue;
						
						if(scaleEstimates[planeIds[i][p1]][cInd].noInliers > scaleEstimates[planeIds[j][p2]][cInd2].noInliers)
						{
							scaleEstimates[planeIds[j][p2]][cInd2].noInliers = scaleEstimates[planeIds[i][p1]][cInd].noInliers;
							scaleEstimates[planeIds[j][p2]][cInd2].scale = 1.0/scaleEstimates[planeIds[i][p1]][cInd].scale;
						}
						else
						{
							scaleEstimates[planeIds[i][p1]][cInd].noInliers = scaleEstimates[planeIds[j][p2]][cInd2].noInliers;
							scaleEstimates[planeIds[i][p1]][cInd].scale = 1.0/scaleEstimates[planeIds[j][p2]][cInd2].scale;
						}
					}
				}
			}
		}
		
		return ;
	}
	
	for(int i=0; i<noImages; i++)
	{
		noPlanes += noRectifiedImages[i];
	}
	
	scaleEstimates.resize(noPlanes);
	
	int curPlane = 0;
	for(int i=0; i<noImages; i++)
	{
		for(int p1=0; p1<noRectifiedImages[i]; p1++)
		{
			char keyFile[1024];
			sprintf(keyFile, "%08d_%d.%s", i, p1, "key");
			string filename = mainDirectory + "/rectification/" + keyFile;
			
			vector<short> refKeys;
			vector<keypt_t> refKeyInfo;
			ifstream keyFin(filename.c_str());
			if(keyFin)
			{
				kpExtractor->ReadKeyFile(filename.c_str(), refKeys, refKeyInfo);
			}
			else
			{
				sift->processUprightImage(rectifiedImages[i][p1], refKeyInfo, refKeys);
				kpExtractor->WriteKeyFile(filename.c_str(), refKeys, refKeyInfo);
			}
			for(int j=0; j<noImages; j++)
			{
				for(int p2=0; p2<noRectifiedImages[j]; p2++)
				{
					if(i==j && p1 == p2)
						continue;
					
					char keyFile[1024];
					sprintf(keyFile, "%08d_%d.%s", j, p2, "key");
					string filename = mainDirectory + "/rectification/" + keyFile;
					
					vector<short> curKeys;
					vector<keypt_t> curKeyInfo;
					ifstream keyFin(filename.c_str());
					if(keyFin)
					{
						kpExtractor->ReadKeyFile(filename.c_str(), curKeys, curKeyInfo);
					}
					else
					{
						sift->processUprightImage(rectifiedImages[j][p2], curKeyInfo, curKeys);
						kpExtractor->WriteKeyFile(filename.c_str(), curKeys, curKeyInfo);
					}
					
					vector<float> potScales;
					double scale = 0.0;
					int noInliers;
					std::vector<KeypointMatch> m = kpMatcher->MatchKeys(curKeyInfo.size(), curKeys, curKeyInfo, refKeyInfo.size(), refKeys, refKeyInfo,
																		threshold, scale, potScales, noInliers);
					scaleEstimate s;
					s.imgIndex1 = i;
					s.imgIndex2 = j;
					s.planeIndex1 = p1;
					s.planeIndex2 = p2;
					s.scale = scale;
					s.noInliers = noInliers;
					scaleEstimates[curPlane].push_back(s);
				}
			}
			curPlane += 1;
		}
	}
	
	int count = 0;
	
	ofstream fout(scaleFilename.c_str(), ios::out);
	fout << noPlanes << endl;
	
	for(int i=0; i<noImages; i++)
	{
		for(int p=0; p<noRectifiedImages[i]; p++)
		{
			sort(scaleEstimates[count].begin(), scaleEstimates[count].end());
			
			fout << i << " " << p << " " << scaleEstimates[count].size() << endl;
			for(int c=0; c<scaleEstimates[count].size(); c++)
			{
				fout << scaleEstimates[count][c].imgIndex2 << " " << scaleEstimates[count][c].planeIndex2 << " " << scaleEstimates[count][c].noInliers << " " << scaleEstimates[count][c].scale << endl;
 			}
			count ++;
		}
	}
	
	fout.close();
}

void BundlerManager::findMatches(int ind, int rectificationInd, Vec2f topLeft, int w, int h)
{
	if(rectifiedImages.size() == 0)
	{
		printf("Rectified image not found.\n");
		return ;
	}
	
	if(rectifiedImages[ind][0] == NULL)
	{
		printf("Rectified image not found.\n");
		return ;
	}
	
	if(ind > imageFileNames.size() -1)
	{
		printf("Image index out of bounds.\n");
		return ;
	}
	
	int noImages = imageFileNames.size();
	
	if(kpMatcher == NULL)
		kpMatcher = new KeypointMatcher();
	if(kpExtractor == NULL)
		kpExtractor = new KeypointExtractor();
	
	SIFT *sift = new SIFT();
	RepetitionFinder rf;
	
	if(!keySpaceAllocated)
	{
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
		keySpaceAllocated = true;
	}
	else
	{
		for(int i=0; i<keys.size(); i++)
			keys[i].clear();
		keys.clear();
		
		for(int i=0; i<keyInfo.size(); i++)
			keyInfo[i].clear();
		keyInfo.clear();
		
		noFeaturePoints.clear();
		
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
	}
	
	for(int i=0; i<matches.size(); i++)
	{
		for(int j=0; j<matches[i].size(); j++)
			matches[i][j].clear();
		matches[i].clear();
	}
	matches.clear();
	matches.resize(noImages);
	for(int i=0; i<noImages; i++)
	{
		matches[i].resize(noImages);
	}
	
	//read previous matches
	readAllGridMatches();
	
	vector<float> scales;
	estimateScales();
	findScaleOrder(ind, rectificationInd, scales);
	
	vector<vector<bool> > matched;
	matched.resize(noImages);
	
	for(int i=0; i<noImages; i++)
		matched[i].resize(noRectifiedImages[i], false);
	
	//matching threshold score
	float thresh = 0.6;
	
	//reference image
	int index = ind;
	Img temp;
	
	Vec2i orgTempSize;
	
	vector<Vec2f> tempMatches;
	vector<float> scores;
	Vec2f gridVectorX, gridVectorY;
	
	vector<vector<Vec2i> > newElementSizes;
	vector<vector<vector<repetitionGrid> > > newRepeatingGrids;
	vector<vector<vector<vector<vector<Vec2f> > > > > newRepetitions;
	
	newRepeatingGrids.resize(noImages);
	newRepetitions.resize(noImages);
	
	newElementSizes.resize(noImages);
	newElementSizes[index].resize(noRectifiedImages[index]);
	newElementSizes[index][rectificationInd] = Vec2i(w,h);
	
	orgTempSize = Vec2i(w,h);
		
	temp.resize(w, h, 0.0);
		
	for(int y=0; y<h; y++)
	{
		for(int x=0; x<w; x++)
		{
			temp.setColor(x, y, (*(rectifiedImages[index][rectificationInd]))(topLeft[0]+x, topLeft[1]+y));
		}
	}
	
	//find matches on reference image
	rf.templateMatchingWithFFTCorrelation(&temp, rectifiedImages[index][rectificationInd], tempMatches, scores, thresh);
	vector<vector<Vec2f> > finalGroups;
	vector<Vec2f> finalTransformations;
	rf.findGridGeneratingVectors(tempMatches, scores, w, h, finalGroups, finalTransformations);
	newRepetitions[index].resize(noRectifiedImages[index]);
	newRepeatingGrids[index].resize(noRectifiedImages[index]);
	
	newRepetitions[index][rectificationInd].resize(finalGroups.size());
	
	matched[ind][rectificationInd] = true;
	
	for(int group=0; group<finalGroups.size(); group++)
	{
		repetitionGrid grid;
		grid.xTrans = Vec2f(finalTransformations[group][0], 0.0);
		grid.yTrans = Vec2f(0.0, finalTransformations[group][1]);
		
		vector<float> groupScores;
		groupScores.resize(finalGroups[group].size(), float(group+1)/(float)(finalGroups.size()));
		rf.segmentIntoRowsAndColumns(finalGroups[group], groupScores, grid, grid.yTrans, grid.xTrans, w, h);
		newRepeatingGrids[index][rectificationInd].push_back(grid);
		
		for(int j=0; j<grid.gridCells.size(); j++)
		{
			Vec2f center = grid.gridCells[j].center;
			
			vector<Vec2f> coord;
			Vec3f pos;
			pos = rectifyingHomographies[index][rectificationInd] * Vec3f(center[0] - w/2, center[1] - h/2, 1.0);
			coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
			pos = rectifyingHomographies[index][rectificationInd] * Vec3f(center[0] + w/2, center[1] - h/2, 1.0);
			coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
			pos = rectifyingHomographies[index][rectificationInd] * Vec3f(center[0] + w/2, center[1] + h/2, 1.0);
			coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
			pos = rectifyingHomographies[index][rectificationInd] * Vec3f(center[0] - w/2, center[1] + h/2, 1.0);
			coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
			
			newRepetitions[index][rectificationInd][group].push_back(coord);
			
			newRepeatingGrids[index][rectificationInd][newRepeatingGrids[index][rectificationInd].size()-1].gridCells[j].size = Vec2i(w,h);
		}
	}
	
	double threshold = 0.05;
	
	int count = 0;
	for(int i=0; i<noImages; i++)
	{
		if(i != ind)
		{
			newElementSizes[i].resize(noRectifiedImages[i]);
			newRepetitions[i].resize(noRectifiedImages[i]);
			newRepeatingGrids[i].resize(noRectifiedImages[i]);
		}
		
		for(int p=0; p<noRectifiedImages[i]; p++)
		{
			if(ind == i && p == rectificationInd)
			{
				count += 1;
				continue;
			}
			
			double scale = scales[count];
			count += 1;
			if(scale == -1)
			{
				newElementSizes[i][p] = Vec2i(0,0);
				continue;
			}
			
			int nW = int(w*scale); int nH = int(h*scale);
			Img nTemp(nW, nH);
			for(int y=0; y<nH; y++)
			{
				for(int x=0; x<nW; x++)
				{
						nTemp.setColor(x, y, temp(x/scale, y/scale));
				}
			}
			newElementSizes[i][p] = Vec2i(nW, nH);
			nTemp.write("/Users/ceylan/Desktop/temp.png");
			vector<Vec2f> nMatches;
			vector<float> nScores;
			rf.templateMatchingWithFFTCorrelation(&nTemp, rectifiedImages[i][p], nMatches, nScores, thresh);
				
			if(nMatches.size() == 0)
				continue;
		
			matched[i][p] = true;
					
			vector<vector<Vec2f> > finalGroups;
			vector<Vec2f> finalTransformations;
			rf.findGridGeneratingVectors(nMatches, nScores, nW, nH, finalGroups, finalTransformations);
				
			newRepetitions[i][p].resize(finalGroups.size());
				
			for(int group=0; group<finalGroups.size(); group++)
			{
				repetitionGrid grid;
				grid.xTrans = Vec2f(finalTransformations[group][0], 0.0);
				grid.yTrans = Vec2f(0.0, finalTransformations[group][1]);
				
				vector<float> groupScores;
				groupScores.resize(finalGroups[group].size(), float(group+1)/(float)(finalGroups.size()));
				rf.segmentIntoRowsAndColumns(finalGroups[group], groupScores, grid, grid.yTrans, grid.xTrans, w, h);
				newRepeatingGrids[i][p].push_back(grid);
					
				for(int j=0; j<grid.gridCells.size(); j++)
				{
					Vec2f center = grid.gridCells[j].center;
					
					vector<Vec2f> coord;
					Vec3f pos;
					pos = rectifyingHomographies[i][p] * Vec3f(center[0] - w/2, center[1] - h/2, 1.0);
					coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
					pos = rectifyingHomographies[i][p] * Vec3f(center[0] + w/2, center[1] - h/2, 1.0);
					coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
					pos = rectifyingHomographies[i][p] * Vec3f(center[0] + w/2, center[1] + h/2, 1.0);
					coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
					pos = rectifyingHomographies[i][p] * Vec3f(center[0] - w/2, center[1] + h/2, 1.0);
					coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
						
					newRepetitions[i][p][group].push_back(coord);
				
					newRepeatingGrids[i][p][newRepeatingGrids[i][p].size()-1].gridCells[j].size = Vec2i(w,h);
				}
			}
		}
	}
	
	/*for(int i=0; i<noImages; i++)
	{
		for(int p=0; p<noRectifiedImages[i]; p++)
		{
			if(matched[i][p])
				continue;
			
			float scale;
			if(true)//saveScale)
			{
				//cluster potential scale values
				int maxClusterSize = -1;
				int maxClusterIndex = -1;
				float lastValue;
			
				for(int loop=0; loop<potentialScales[i][p].size(); loop++)
				{
					lastValue = potentialScales[i][p][loop];
					int count = 0;
					for(int j=0; j<potentialScales[i][p].size(); j++)
					{
						if(abs(potentialScales[i][p][j] - lastValue) <= threshold)
						{
							count++;
						}
					}
					if(maxClusterSize == -1)
					{
						maxClusterSize = count;
						maxClusterIndex = loop;
					}
					else if(count > maxClusterSize)
					{
						maxClusterSize = count;
						maxClusterIndex = loop;
					}
				}
				if(maxClusterIndex != -1 && maxClusterSize >= 40)
				{
					lastValue = potentialScales[i][p][maxClusterIndex];
					vector<float> cluster;
					for(int j=0; j<potentialScales[i][p].size(); j++)
					{
						if(potentialScales[i][p][j] - lastValue <= threshold)
						{
							cluster.push_back(potentialScales[i][p][j]);
						}
					}
					sort(cluster.begin(), cluster.end());
					scale = cluster[cluster.size()/2];
					masks[i][p] = scale;
				}
				else
				{
					scale = masks[i][p];
					if(scale == -1)
					{
						newElementSizes[i][p] = Vec2i(-1, -1);
						continue;
					}
				}
			
				int nW = int(w*scale); int nH = int(h*scale);
				Img nTemp(nW, nH);
				for(int y=0; y<nH; y++)
				{
					for(int x=0; x<nW; x++)
					{
						nTemp.setColor(x, y, temp(x/scale, y/scale));
					}
				}
				newElementSizes[i][p] = Vec2i(nW, nH);
				
				vector<Vec2f> nMatches;
				vector<float> nScores;
			
				rf.templateMatchingWithFFTCorrelation(&nTemp, rectifiedImages[i][p], nMatches, nScores, thresh);
				if(nMatches.size() == 0)
					continue;
				
				vector<vector<Vec2f> > finalGroups;
				vector<Vec2f> finalTransformations;
				rf.findGridGeneratingVectors(nMatches, nScores, nW, nH, finalGroups, finalTransformations);
				
				newRepetitions[i][p].resize(finalGroups.size());
				
				for(int group=0; group<finalGroups.size(); group++)
				{
					repetitionGrid grid;
					grid.xTrans = Vec2f(finalTransformations[group][0], 0.0);
					grid.yTrans = Vec2f(0.0, finalTransformations[group][1]);
					
					vector<float> groupScores;
					groupScores.resize(finalGroups[group].size(), float(group+1)/(float)(finalGroups.size()));
					rf.segmentIntoRowsAndColumns(finalGroups[group], groupScores, grid, grid.yTrans, grid.xTrans, w, h);
					newRepeatingGrids[i][p].push_back(grid);
					
					for(int j=0; j<grid.gridCells.size(); j++)
					{
						Vec2f center = grid.gridCells[j].center;
						
						vector<Vec2f> coord;
						Vec3f pos;
						pos = rectifyingHomographies[i][p] * Vec3f(center[0] - w/2, center[1] - h/2, 1.0);
						coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
						pos = rectifyingHomographies[i][p] * Vec3f(center[0] + w/2, center[1] - h/2, 1.0);
						coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
						pos = rectifyingHomographies[i][p] * Vec3f(center[0] + w/2, center[1] + h/2, 1.0);
						coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
						pos = rectifyingHomographies[i][p] * Vec3f(center[0] - w/2, center[1] + h/2, 1.0);
						coord.push_back(Vec2f(pos[0]/pos[2], pos[1]/pos[2]));
						
						newRepetitions[i][p][group].push_back(coord);
						
						newRepeatingGrids[i][p][newRepeatingGrids[i][p].size()-1].gridCells[j].size = Vec2i(w,h);
					}
				}
			}
		}
	}*/
	
	elementSizes.push_back(newElementSizes);
	repetitions.push_back(newRepetitions);
	repeatingGrids.push_back(newRepeatingGrids);
	
	/*char matchFile[1024];
	sprintf(matchFile, "%08d.%s", ind, "txt");
	string filename = mainDirectory + "/matches/" + matchFile;
	ofstream fout(filename.c_str(), ios::out);
	fout << "1" << endl;
	for(int p=0; p<noRectifiedImages[ind]; p++)
	{
		int noMatches = 0;
		for(int g=0; g<newRepetitions[ind][p].size(); g++)
			noMatches += newRepetitions[ind][p][g].size();
		fout << noMatches << " " << newElementSizes[ind][p][0] << " " << newElementSizes[ind][p][1] << endl;
		for(int g=0; g<newRepeatingGrids[ind][p].size(); g++)
		{
			for(int c=0; c<newRepeatingGrids[ind][p][g].gridCells.size(); c++)
			{
				fout << newRepeatingGrids[ind][p][g].gridCells[c].center[0] << " " << newRepeatingGrids[ind][p][g].gridCells[c].center[1] << " "
					<< newRepeatingGrids[ind][p][g].gridCells[c].descriptor << endl;
				
			}
		}
	}
	fout.close();*/
	
	int noTemplates = elementSizes.size();
	for(int i=0; i<noImages; i++)
	{
		char matchFile[1024];
		sprintf(matchFile, "%08d.%s", i, "txt");
		string filename = mainDirectory + "/matches/" + matchFile;
		ofstream fout(filename.c_str(), ios::out);
		fout << noTemplates << endl;
		for(int t=0; t<noTemplates; t++)
		{
			for(int p=0; p<noRectifiedImages[i]; p++)
			{
				int noMatches = 0;
				for(int g=0; g<repetitions[t][i][p].size(); g++)
					noMatches += repetitions[t][i][p][g].size();
				fout << noMatches << " " << elementSizes[t][i][p][0] << " " << elementSizes[t][i][p][1] << endl;
				for(int g=0; g<repeatingGrids[t][i][p].size(); g++)
				{
					for(int c=0; c<repeatingGrids[t][i][p][g].gridCells.size(); c++)
					{
						fout << repeatingGrids[t][i][p][g].gridCells[c].center[0] << " " << repeatingGrids[t][i][p][g].gridCells[c].center[1] << " "
							<< repeatingGrids[t][i][p][g].gridCells[c].descriptor << endl;
					}
				}
			}
		}
		fout.close();
	}
}


std::vector<std::vector<Vec2f> > BundlerManager::getRepetitions(int index)
{
	vector<vector<Vec2f> > tmpRepetitions;
	
	for(int i=0; i<repetitions.size(); i++)
	{
		for(int j=0; j<noRectifiedImages[index]; j++)
		{
			for(int g=0; g<repetitions[i][index][j].size(); g++)
			{
				for(int k=0; k<repetitions[i][index][j][g].size(); k++)
					tmpRepetitions.push_back(repetitions[i][index][j][g][k]);
			}
		}
	}
	return tmpRepetitions;
}

std::vector<std::vector<Vec2f> > BundlerManager::getRefinedRepetitions(int index)
{
	vector<vector<Vec2f> > tmpRepetitions;
	for(int i=0; i<refinedRepetitions.size(); i++)
	{
		for(int j=0; j<refinedRepetitions[i][index].size(); j++)
		{
			for(int g=0; g<refinedRepetitions[i][index][j].size(); g++)
			{
				tmpRepetitions.push_back(refinedRepetitions[i][index][j][g]);
			}
		}
	}
	return tmpRepetitions;
}

std::vector<Color> BundlerManager::getRepetitionColors(int index)
{
	vector<Color> scores;
	for(int t = 0; t<repeatingGrids.size(); t++)
	{
		for(int p=0; p<repeatingGrids[t][index].size(); p++)
		{
			for(int g=0; g<repeatingGrids[t][index][p].size(); g++)
			{
				for(int i=0; i<repeatingGrids[t][index][p][g].gridCells.size(); i++)
				{
					Vec3uc col = MathUtils::generateColorFromValue(repeatingGrids[t][index][p][g].gridCells[i].descriptor, 0.0, 1.0);
					scores.push_back(Color(col[0]/255.0, col[1]/255.0, col[2]/255.0));
				}
			}
		}
	}
	return scores;
}

void BundlerManager::createRepetitionMask(int index)
{
	/*char imgFile[1024];
	sprintf(imgFile, "%08d.%s", index, "png");
	string filename = mainDirectory + "/masks/" + imgFile;
	
	Img currentImage;
	currentImage.read(imageFileNames[index]);
	Img maskImage;
	maskImage.resize(currentImage.width(), currentImage.height(), Color(1.0, 1.0, 1.0));
	
	for(int i=0; i<repeatingGrids[index].gridCells.size(); i++)
	{
		if(repeatingGrids[index].gridCells[i].label == -1)
			continue;
		
		Vec2f center = repeatingGrids[index].gridCells[i].center;
		Vec2i size = repeatingGrids[index].gridCells[i].size;
		int offsetX = max(0, (int)(repeatingGrids[index].xTrans[0]) - size[0])/2;
		int offsetY = max(0, (int)(repeatingGrids[index].yTrans[1]) - size[1])/2;
		size[0] = size[0] + offsetX;
		size[1] = size[1] + offsetY;
		
		for(int x=center[0]-size[0]/2; x<center[0]+size[0]/2; x++)
		{
			for(int y=center[1]-size[1]/2; y<center[1]+size[1]/2; y++)
			{
				Vec3f pos = rectifyingHomographies[index] * Vec3f(x,y,1.0);
				if(pos[2] != 0.0)
				{
					pos[0] /= pos[2]; pos[1] /= pos[2];
					maskImage.setColor(pos[0], pos[1], repetitionMaskColors[repeatingGrids[index].gridCells[i].label]);
				}
			}
		}
	}
	
	maskImage.write(filename.c_str());*/
}

void BundlerManager::ComputeAllFeaturePoints()
{	
	SIFT *sift = new SIFT();
	
	if(kpExtractor == NULL)
		kpExtractor = new KeypointExtractor();

	if(!keySpaceAllocated)
	{
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
		keySpaceAllocated = true;
	}
	else
	{
		for(int i=0; i<keys.size(); i++)
			keys[i].clear();
		keys.clear();
		
		for(int i=0; i<keyInfo.size(); i++)
			keyInfo[i].clear();
		keyInfo.clear();
		
		noFeaturePoints.clear();
		
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
	}
	
	for(int i=0; i<imageFileNames.size(); i++)
	{
		string filename = imageFileNames[i];
		string keyFilename = filename.substr(0, filename.find(".",0));
		keyFilename += ".key";
		
        printf("keyFile:%s\n", keyFilename.c_str());
		ifstream fin(keyFilename.c_str());
		if(fin)
			ReadFeaturePoints(i);
		else
		{
			Img im;
			im.read(imageFileNames[i].c_str());
			sift->processImage(&im, keyInfo[i], keys[i]);
			noFeaturePoints[i] = keyInfo[i].size();
		}
		
	}
}

void BundlerManager::ReadFeaturePoints(int imageIndex)
{
	string filename = imageFileNames[imageIndex];
	string keyFilename = filename.substr(0, filename.find(".",0));
	keyFilename += ".key";
	
	if(kpExtractor == NULL)
		kpExtractor = new KeypointExtractor();
	
	int n = kpExtractor->ReadKeyFile(keyFilename.c_str(), keys[imageIndex], keyInfo[imageIndex]);
	noFeaturePoints[imageIndex] = n;
	
	printf("Image %d has %d keypoints.\n", imageIndex, n);
}

void BundlerManager::ReadFeaturePointsWithRepetition(int imageIndex)
{
	char keyFile[1024];
	sprintf(keyFile, "%08d.%s", imageIndex, "key");
	string keyFilename = mainDirectory + "/matches/" + keyFile;
	
	if(kpExtractor == NULL)
		kpExtractor = new KeypointExtractor();
	
	int n = kpExtractor->ReadKeyFile(keyFilename.c_str(), keys[imageIndex], keyInfo[imageIndex]);
	noFeaturePoints[imageIndex] = n;
	
	printf("Image %d has %d keypoints.\n", imageIndex, n);
}

void BundlerManager::ReadAllFeaturePoints()
{
	if(kpExtractor == NULL)
		kpExtractor = new KeypointExtractor();
	
	if(!keySpaceAllocated)
	{
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
		keySpaceAllocated = true;
	}
	else
	{
		for(int i=0; i<keys.size(); i++)
			keys[i].clear();
		keys.clear();
		
		for(int i=0; i<keyInfo.size(); i++)
			keyInfo[i].clear();
		keyInfo.clear();
		
		noFeaturePoints.clear();
		
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
	}
	
	for(int i=0; i<imageFileNames.size(); i++)
	{
		ReadFeaturePoints(i);
	}
	
	//formMatchesFromTracks();
}

void BundlerManager::ReadAllFeaturePointsWithRepetition()
{
	if(kpExtractor == NULL)
		kpExtractor = new KeypointExtractor();
	
	if(!keySpaceAllocated)
	{
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
		keySpaceAllocated = true;
	}
	else
	{
		for(int i=0; i<keys.size(); i++)
			keys[i].clear();
		keys.clear();
		
		for(int i=0; i<keyInfo.size(); i++)
			keyInfo[i].clear();
		keyInfo.clear();
		
		noFeaturePoints.clear();
		
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
	}
	
	for(int i=0; i<imageFileNames.size(); i++)
	{
		ReadFeaturePointsWithRepetition(i);
	}

}

void BundlerManager::findUniqueFeaturePoints()
{
	int noImages = imageFileNames.size();
	for(int i=0; i<noImages; i++)
	{
		vector<KeypointMatch> mat = kpMatcher->findUniqueKeys(noFeaturePoints[i], keys[i]);
	 
		vector<keypt_t> keyInfoNew;
		vector<short> keysNew; 
		for(int m=0; m<mat.size(); m++)
		{
			int idx = mat[m].m_idx1;
			keyInfoNew.push_back(keyInfo[i][idx]);
			for(int d=0; d<128; d++)
			{
				keysNew.push_back(keys[i][idx*128+d]);
			}
		}
		printf("keys before:%d, keys after:%d\n", keyInfo[i].size(), keyInfoNew.size());
		keyInfo[i].clear();
		keys[i].clear();
		for(int m=0; m<keyInfoNew.size(); m++)
		{
			keyInfo[i].push_back(keyInfoNew[m]);
			for(int d=0; d<128; d++)
			{
				keys[i].push_back(keysNew[m*128+d]);
			}
		}
		noFeaturePoints[i] = keyInfo[i].size();
	}
	
}

void BundlerManager::getFeaturePointInfo(int imageIndex, keypt_t* &featurePointInfo, int &noFeatures)
{
	if(imageIndex >	noFeaturePoints.size()-1 || !keySpaceAllocated)
	{
		printf("Not valid image index.\n");
		return;
	}
	
	/*noFeatures = keyPairs[imageIndex].size();
	featurePointInfo = new keypt_t[keyPairs[imageIndex].size()*2];
	for(int i=0; i<keyPairs[imageIndex].size(); i++)
	{
		featurePointInfo[i*2] = keyPairs[imageIndex][i].keyInfo1;
		featurePointInfo[i*2+1] = keyPairs[imageIndex][i].keyInfo2;
	}*/
	
	noFeatures = noFeaturePoints[imageIndex];
	if(noFeatures > 0)
		featurePointInfo = &(keyInfo[imageIndex][0]);
}

void BundlerManager::addFeaturePoint(int imageIndex, Vec2f pt)
{
	if(!keySpaceAllocated)
	{
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
		keySpaceAllocated = true;
	}
	/*else
	{
		for(int i=0; i<keys.size(); i++)
			keys[i].clear();
		keys.clear();
		
		for(int i=0; i<keyInfo.size(); i++)
			keyInfo[i].clear();
		keyInfo.clear();
		
		noFeaturePoints.clear();
		
		keys.resize(imageFileNames.size());
		keyInfo.resize(imageFileNames.size());
		noFeaturePoints.resize(imageFileNames.size(), 0);
	}*/
	
	for(int i=0; i<128; i++)
		keys[imageIndex].push_back(0);
	
	keypt_t kp;
	kp.x = pt[0];
	kp.y = pt[1];
	kp.scale = 0.0;
	kp.orient = 0.0;
	keyInfo[imageIndex].push_back(kp);
	
	noFeaturePoints[imageIndex] += 1;
}

void BundlerManager::saveFeaturePoints()
{
	if(keySpaceAllocated)
	{
		int numImages = imageFileNames.size();
		for(int i=0; i<numImages; i++)
		{
			string fname = imageFileNames[i];
			fname = fname.substr(0, fname.find("."));
			fname += ".key";
			FILE *f = fopen(fname.c_str(), "w");
			fprintf(f, "%d %d\n", noFeaturePoints[i], 128);
			for(int j=0; j<noFeaturePoints[i]; j++)
			{
				fprintf(f, "%f %f %f %f\n", keyInfo[i][j].y, keyInfo[i][j].x, keyInfo[i][j].scale, keyInfo[i][j].orient);
				for(int l=0;l<7; l++)
				{
					if(l<6)
					{
						fprintf(f, "%hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd\n", 
								keys[i][j*128+l*20+0], keys[i][j*128+l*20+1], keys[i][j*128+l*20+2], keys[i][j*128+l*20+3], keys[i][j*128+l*20+4],
								keys[i][j*128+l*20+5], keys[i][j*128+l*20+6], keys[i][j*128+l*20+7], keys[i][j*128+l*20+8], keys[i][j*128+l*20+9],
								keys[i][j*128+l*20+10], keys[i][j*128+l*20+11], keys[i][j*128+l*20+12], keys[i][j*128+l*20+13], keys[i][j*128+l*20+14],
								keys[i][j*128+l*20+15], keys[i][j*128+l*20+16], keys[i][j*128+l*20+17], keys[i][j*128+l*20+18], keys[i][j*128+l*20+19]);
						
					}
					else
					{
						fprintf(f,"%hd %hd %hd %hd %hd %hd %hd %hd\n",
								keys[i][j*128+l*20+0], keys[i][j*128+l*20+1], keys[i][j*128+l*20+2], keys[i][j*128+l*20+3],
								keys[i][j*128+l*20+4], keys[i][j*128+l*20+5], keys[i][j*128+l*20+6], keys[i][j*128+l*20+7]);
					}
				}
			}
			fclose(f);
		}
	}
}

void BundlerManager::buildManualCorrespondences()
{
	int numImages = imageFileNames.size();
	for(int i=0; matches.size(); i++)
	{
		for(int j=0; j<matches[i].size(); j++)
			matches[i][j].clear();
		matches[i].clear();
	}
	matches.clear();
	matches.resize(numImages);
	
	for(int i=0; i<numImages; i++)
	{
		int noFeatures = noFeaturePoints[i];
		
		for (int j = 0; j < i; j++) 
		{
			vector<KeypointMatch> m;
			
			for(int k=0; k<noFeatures; k++)
			{
				if(k == keys[j].size())
					break;
				
				Vec2f p1 = Vec2f(keyInfo[i][k].x, keyInfo[i][k].y);
				Vec2f p2 = Vec2f(keyInfo[j][k].x, keyInfo[j][k].y);
				
				if(p1 != Vec2f(0.0, 0.0) && p2 != Vec2f(0.0, 0.0))
				{
					KeypointMatch kp(k, k);
					kp.validMatch = true;
					m.push_back(kp);
				}
			}
			matches[j].push_back(m);
		}
	}
}

void BundlerManager::matchFeaturePoints(bool useMask)
{
	if(kpMatcher == NULL)
		kpMatcher = new KeypointMatcher();
	
	int numImages = imageFileNames.size();
	double ratio = 0.6;
	
		
	matches.clear();
	matches.resize(numImages);
	
	if(maskFileNames.size() != numImages)
	{
		printf("Mask images not specified.\n");
		useMask = false;
	}
	
	for (int i = 0; i < numImages; i++) 
	{
		if (noFeaturePoints[i] == 0)
			continue;
		
		for (int j = 0; j < i; j++) 
		{
			std::vector<KeypointMatch> m;
		
			if (noFeaturePoints[j] != 0)
			{
				printf("Matching image %d to image %d\n", j, i);
			
				// Compute likely matches between two sets of keypoints
				if(useMask)
				{
					Img img1, img2;
					img1.read(maskFileNames[j]);
					img2.read(maskFileNames[i]);
					
					std::vector<KeypointMatch> m = kpMatcher->MatchKeys(noFeaturePoints[j], keys[j], keyInfo[j], noFeaturePoints[i], keys[i], keyInfo[i], ratio, img1, img2);
					
					std::vector<KeypointMatch> finalMatches;
					for(int k=0; k<m.size(); k++)
					{
						keypt_t kp1 = keyInfo[j][m[k].m_idx1];
						keypt_t kp2 = keyInfo[i][m[k].m_idx2];
						
						//if(img1(kp1.x, kp1.y) == img2(kp2.x, kp2.y) /*&& img1(kp1.x, kp1.y) != Color(1.0, 1.0, 1.0) && img2(kp2.x, kp2.y) != Color(1.0, 1.0, 1.0)*/)
							finalMatches.push_back(m[k]);
					}
					matches[j].push_back(finalMatches);
				}
				else
				{
					m = kpMatcher->MatchKeys(noFeaturePoints[j], keys[j], noFeaturePoints[i], keys[i], ratio);
					matches[j].push_back(m);
				}
			}
			else
			{
				matches[j].push_back(std::vector<KeypointMatch>());
			}
		}
	}
	
	computeTracks();
}

void BundlerManager::writeMatches()
{
	string dir = mainDirectory;
	string filename = dir + "/matches_init.txt";
	printf("file:%s\n", filename.c_str());
	
	if(kpMatcher == NULL)
		kpMatcher = new KeypointMatcher();
	
	kpMatcher->writeMatches(filename.c_str(), matches);
	
}

void BundlerManager::readMatches()
{
	string dir = mainDirectory;
	string filename = dir + "/matches_init.txt";
	printf("file:%s\n", filename.c_str());
	
	if(kpMatcher == NULL)
		kpMatcher = new KeypointMatcher();
	
	matches.clear();
	matches.resize(imageFileNames.size());
	fundMatrices.resize(imageFileNames.size());
	for(int i=0; i<imageFileNames.size(); i++)
	{
		matches[i].resize(imageFileNames.size()-i-1);
		fundMatrices[i].resize(imageFileNames.size()-i-1);
	}
	
	kpMatcher->readMatches(filename.c_str(), matches);
	
	EpipolarGeometry eg;
	
	string fundFile = dir + "/fundMatrix.txt";
	ifstream fin(fundFile.c_str(), ios::in);
	if(fin)
	{
		int noMatrices;
		fin >> noMatrices;
		
		for(int i=0; i<imageFileNames.size(); i++)
		{
			for(int j=i+1; j<imageFileNames.size(); j++)
			{
				Matrix3f F;
				F.setIdentity();
				
				fin >> F[0][0] >> F[1][0] >> F[2][0];
				fin >> F[0][1] >> F[1][1] >> F[2][1];
				fin >> F[0][2] >> F[1][2] >> F[2][2];
				
				printf("Fund matrix for image %d and %d:\n", i, j);
				printf("%f %f %f\n", F[0][0], F[1][0], F[2][0]);
				printf("%f %f %f\n", F[0][1], F[1][1], F[2][1]);
				printf("%f %f %f\n", F[0][2], F[1][2], F[2][2]);
				
				fundMatrices[i][j-i-1] = F;
			}
		}
	}
	else
	{
		ofstream fout(fundFile.c_str(), ios::out);
		fout << (imageFileNames.size()*(imageFileNames.size()-1)/2) << endl;
	
		for(int i=0; i<imageFileNames.size(); i++)
		{
			for(int j=i+1; j<imageFileNames.size(); j++)
			{
				printf("matching image %d to image %d\n", i, j);
				Matrix3f F;
				F.setIdentity();
				if(matches[i][j-i-1].size() >= 8)
					eg.estimateFMatrix(keyInfo[i], keyInfo[j], matches[i][j-i-1], 1024, 20.0, F);
				fout << F[0][0] << " " << F[1][0] << " " << F[2][0] << endl;
				fout << F[0][1] << " " << F[1][1] << " " << F[2][1] << endl;
				fout << F[0][2] << " " << F[1][2] << " " << F[2][2] << endl;
				fundMatrices[i][j-i-1] = F;
			}
		}
		fout.close();
	}
	
	//ucl1:cams 24,26 matches:10,11,19,21
	//ucl5:cams 0,1 matches:16,18,22,23
	//ucl5lowRes cams 0,1 matches:33,36,42,44
	//ucl9:cams 3 4 matches:14,15,17,18
	//museum:cam 0,1 matches:12,13,26,28
	//pub: cam 0,1 matches:48,59,34,42
	//ucl3:cam 0,1 matches:48,50,54,58
    //lausanne1: cam 0,1 matches 30,33,48,49
	vector<Vec2f> imageCoordinates1;
	vector<Vec2f> imageCoordinates2;
	int cam1 = 0;
	int cam2 = 1;
	imageCoordinates1.push_back(Vec2f(keyInfo[cam1][matches[cam1][cam2-cam1-1][16].m_idx1].x, keyInfo[cam1][matches[cam1][cam2-cam1-1][16].m_idx1].y));
	imageCoordinates1.push_back(Vec2f(keyInfo[cam1][matches[cam1][cam2-cam1-1][18].m_idx1].x, keyInfo[cam1][matches[cam1][cam2-cam1-1][18].m_idx1].y));
	imageCoordinates1.push_back(Vec2f(keyInfo[cam1][matches[cam1][cam2-cam1-1][22].m_idx1].x, keyInfo[cam1][matches[cam1][cam2-cam1-1][22].m_idx1].y));
	imageCoordinates1.push_back(Vec2f(keyInfo[cam1][matches[cam1][cam2-cam1-1][23].m_idx1].x, keyInfo[cam1][matches[cam1][cam2-cam1-1][23].m_idx1].y));
	
	imageCoordinates2.push_back(Vec2f(keyInfo[cam2][matches[cam1][cam2-cam1-1][16].m_idx2].x, keyInfo[cam2][matches[cam1][cam2-cam1-1][16].m_idx2].y));
	imageCoordinates2.push_back(Vec2f(keyInfo[cam2][matches[cam1][cam2-cam1-1][18].m_idx2].x, keyInfo[cam2][matches[cam1][cam2-cam1-1][18].m_idx2].y));
	imageCoordinates2.push_back(Vec2f(keyInfo[cam2][matches[cam1][cam2-cam1-1][22].m_idx2].x, keyInfo[cam2][matches[cam1][cam2-cam1-1][22].m_idx2].y));
	imageCoordinates2.push_back(Vec2f(keyInfo[cam2][matches[cam1][cam2-cam1-1][23].m_idx2].x, keyInfo[cam2][matches[cam1][cam2-cam1-1][23].m_idx2].y));
	
	//bringToSameCoordinateFrame(imageCoordinates1, imageCoordinates2, cam1, cam2);
     
}

int BundlerManager::getNoMatches(int imageIndex1, int imageIndex2)
{
	if(imageIndex1 < imageIndex2)
	{
		return matches[imageIndex1][imageIndex2-imageIndex1-1].size();
	}
	else
	{
		return matches[imageIndex2][imageIndex1-imageIndex2-1].size();
	}
}

void BundlerManager::getMatchingPoints(int imageIndex1, int imageIndex2, int matchIndex, Vec2f &pt1, Vec2f &pt2, Vec2f &epiStart, Vec2f &epiEnd, bool &valid)
{
	KeypointMatch kp;
	
	epiStart = Vec2f(0.0, 0.0);
	epiEnd = Vec2f(0.0, 0.0);
	
	EpipolarGeometry eg;
	
	if(imageIndex1 < imageIndex2)
	{
		if(matchIndex < matches[imageIndex1][imageIndex2-imageIndex1-1].size())
		{
			kp = matches[imageIndex1][imageIndex2-imageIndex1-1][matchIndex];
			pt1[0] = keyInfo[imageIndex1][kp.m_idx1].x;
			pt1[1] = keyInfo[imageIndex1][kp.m_idx1].y;
			pt2[0] = keyInfo[imageIndex2][kp.m_idx2].x;
			pt2[1] = keyInfo[imageIndex2][kp.m_idx2].y;
			
			if(transformations.size() > 0)
			{
				printf("F:\n");
				Matrix3f F = (transformations[imageIndex1][imageIndex2]).fMatrix;
				for(int i=0; i<3; i++)
				{
					for(int j=0; j<3; j++)
					{
						printf("%f ", F[j][i]);
					}
					printf("\n");
				}
				Vec3f l = (transformations[imageIndex1][imageIndex2]).fMatrix * Vec3f(pt1[0], pt1[1], 1.0);
				GeometryUtils::computeBoundariesOfTheLine(l, 2000, 1500, epiStart, epiEnd);
			}
			
			/*float dist = eg.fmatrixComputeResidual(fundMatrices[imageIndex1][imageIndex2-imageIndex1-1], expand2To3(pt2), expand2To3(pt1));
			if(dist < 50.0)
				valid = true;
			else
				valid = false;
			
			Vec3f l = fundMatrices[imageIndex1][imageIndex2-imageIndex1-1] * Vec3f(pt1[0], pt1[1], 1.0);
			GeometryUtils::computeBoundariesOfTheLine(l, 2160, 2880, epiStart, epiEnd);
			valid = true;*/
		}
	}
	else
	{
		if(matchIndex < matches[imageIndex2][imageIndex1-imageIndex2-1].size())
		{
			kp = matches[imageIndex2][imageIndex1-imageIndex2-1][matchIndex];
			pt2[0] = keyInfo[imageIndex2][kp.m_idx1].x;
			pt2[1] = keyInfo[imageIndex2][kp.m_idx1].y;
			pt1[0] = keyInfo[imageIndex1][kp.m_idx2].x;
			pt1[1] = keyInfo[imageIndex1][kp.m_idx2].y;
			
			if(transformations.size() > 0)
			{
				printf("F:\n");
				Matrix3f F = (transformations[imageIndex1][imageIndex2]).fMatrix;
				for(int i=0; i<3; i++)
				{
					for(int j=0; j<3; j++)
					{
						printf("%f ", F[j][i]);
					}
					printf("\n");
				}
				
				Vec3f l = (transformations[imageIndex1][imageIndex2]).fMatrix * Vec3f(pt1[0], pt1[1], 1.0);
				GeometryUtils::computeBoundariesOfTheLine(l, 2160, 2880, epiStart, epiEnd);
			}
			
			/*float dist = eg.fmatrixComputeResidual(fundMatrices[imageIndex2][imageIndex1-imageIndex2-1], expand2To3(pt1), expand2To3(pt2));
			if(dist < 50.0)
				valid = true;
			else
				valid = false;
			
			Vec3f l = fundMatrices[imageIndex2][imageIndex1-imageIndex2-1].transpose() * Vec3f(pt1[0], pt1[1], 1.0);
			GeometryUtils::computeBoundariesOfTheLine(l, 2160, 2880, epiStart, epiEnd);
			valid = true;*/
		}
	}
	
}

void BundlerManager::getAllMatchingPoints(int imageIndex1, int imageIndex2, vector<Vec2f> &pt1, vector<Vec2f> &pt2, vector<bool> &valid)
{
	KeypointMatch kp;
	
	EpipolarGeometry eg;
	
		
	if(imageIndex1 < imageIndex2)
	{
		for(int i=0; i<matches[imageIndex1][imageIndex2-imageIndex1-1].size(); i++)
		{
			kp = matches[imageIndex1][imageIndex2-imageIndex1-1][i];
			pt1.push_back(Vec2f(keyInfo[imageIndex1][kp.m_idx1].x, keyInfo[imageIndex1][kp.m_idx1].y));
			pt2.push_back(Vec2f(keyInfo[imageIndex2][kp.m_idx2].x, keyInfo[imageIndex2][kp.m_idx2].y));
			
			float dist = eg.fmatrixComputeResidual(fundMatrices[imageIndex1][imageIndex2-imageIndex1-1], expand2To3(pt2[i]), expand2To3(pt1[i]));
			if(dist < 50.0)
				valid.push_back(true);
			else
				valid.push_back(false);
		}
	}
	else
	{
		for(int i=0; i<matches[imageIndex2][imageIndex1-imageIndex2-1].size(); i++)
		{
			kp = matches[imageIndex2][imageIndex1-imageIndex2-1][i];
			pt2.push_back(Vec2f(keyInfo[imageIndex2][kp.m_idx1].x, keyInfo[imageIndex2][kp.m_idx1].y));
			pt1.push_back(Vec2f(keyInfo[imageIndex1][kp.m_idx2].x, keyInfo[imageIndex1][kp.m_idx2].y));
			
			float dist = eg.fmatrixComputeResidual(fundMatrices[imageIndex2][imageIndex1-imageIndex2-1], expand2To3(pt1[i]), expand2To3(pt2[i]));
			if(dist < 50.0)
				valid.push_back(true);
			else
				valid.push_back(false);
		}
	}
	
}

float* BundlerManager::getProjectedPoints(int &noPoints, int imageIndex)
{
	if(pointProjections.size() > imageIndex)
	{
		noPoints = pointProjections[imageIndex].size();
		return &(pointProjections[imageIndex][0][0]);
	}
	else
	{
		noPoints = 0;
		return NULL;
	}
}

Color* BundlerManager::getProjectedPointColors(int &noPoints, int imageIndex)
{
	if(pointProjectionColors.size() > imageIndex)
	{
		noPoints = pointProjectionColors[imageIndex].size();
		return &(pointProjectionColors[imageIndex][0]);
	}
	else
	{
		noPoints = 0;
		return NULL;
	}
}

void BundlerManager::setMatchValidity(int imageIndex1, int imageIndex2, int matchIndex, bool validity)
{
	KeypointMatch kp;
	if(imageIndex1 < imageIndex2)
	{
		matches[imageIndex1][imageIndex2-imageIndex1-1][matchIndex].validMatch = validity;
	}
	else
	{
		matches[imageIndex2][imageIndex1-imageIndex2-1][matchIndex].validMatch = validity;
	}	
}

/*void BundlerManager::writeMatchValidity()
{
	string filename = mainDirectory + "/matchValidity.txt";
	ofstream fout(mainDirectory);
	for(int i=0; i<imageFileNames.size(); i++)
	{
		for(int j=i+1; j<imageFileNames.size(); j++)
		{
			
		}
	}
}*/

void BundlerManager::blendImagesWithMatches(int imageInd1, int imageInd2)
{
	Img imgA, imgB;
	imgA.read(imageFileNames[imageInd1].c_str());
	imgB.read(imageFileNames[imageInd2].c_str());
	
	int width = imgA.width()/4;
	int height = imgA.height()/4;
	int totalStep = 20;
	
	vector<Vec2f> matchInitPos;
	vector<Vec2f> matchDisp;
	
	if(imageInd1 < imageInd2)
	{
		for(int i=0; i<matches[imageInd1][imageInd2-imageInd1-1].size(); i++)
		{
			KeypointMatch m = matches[imageInd1][imageInd2-imageInd1-1][i];
			Vec2f posA = Vec2f(keyInfo[imageInd1][m.m_idx1].x, keyInfo[imageInd1][m.m_idx1].y);
			Vec2f posB = Vec2f(keyInfo[imageInd2][m.m_idx2].x, keyInfo[imageInd2][m.m_idx2].y);
			matchInitPos.push_back(posA);
			matchDisp.push_back((posB-posA)/totalStep);
		}
	}
	else
	{
		for(int i=0; i<matches[imageInd2][imageInd1-imageInd2-1].size(); i++)
		{
			KeypointMatch m = matches[imageInd2][imageInd1-imageInd2-1][i];
			Vec2f posA = Vec2f(keyInfo[imageInd1][m.m_idx2].x, keyInfo[imageInd1][m.m_idx2].y);
			Vec2f posB = Vec2f(keyInfo[imageInd2][m.m_idx1].x, keyInfo[imageInd2][m.m_idx1].y);
			matchInitPos.push_back(posA);
			matchDisp.push_back((posB-posA)/totalStep);
		}
	}
	
	for(int i=0; i<totalStep+1; i++)
	{
		Img tmp(width, height);
#pragma omp parallel for
		for(int x=0; x<width; x++)
		{
			for(int y=0; y<height; y++)
			{
				float alpha = float(totalStep-i)/(float)(totalStep);
				float beta = (float)(i)/(float)(totalStep);
				tmp.setColor(x, y, imgA(x*4,y*4)*alpha + imgB(x*4,y*4)*beta);
			}
		}
		for(int m=0; m<matchInitPos.size(); m++)
		{
			printf("%f %f\n", matchInitPos[m][0]+i*matchDisp[m][0], matchInitPos[m][1]+i*matchDisp[m][1]);
			int posX = floor((matchInitPos[m][0]+i*matchDisp[m][0])/4.0+0.5);
			int posY = floor((matchInitPos[m][1]+i*matchDisp[m][1])/4.0+0.5);
			for(int a=-1; a<2; a++)
			{
				for(int b=-1; b<2; b++)
				{
					if(posX+a>=0 && posX+a<width && posY+b>=0 && posY+b<height)
					{
						tmp.setColor(posX+a, posY+b, Color(1.0,0.0,0.0));
					}
				}
			}
			//tmp.setColor((matchInitPos[m][0]+i*matchDisp[m][0])/4.0, (matchInitPos[m][1]+i*matchDisp[m][1])/4.0, Color(1.0, 0.0, 0.0));
		}
		
		stringstream ss;
		ss << mainDirectory << "/blended_" << i << ".png";
		string filename;
		ss >> filename;
		tmp.write(filename.c_str());
	}
}

void BundlerManager::registerImages()
{
	int numImages = imageFileNames.size();
	GraphWrapper gw(numImages);
	vector<vector<vector<imageTransformation> > > imageAlignments;
	
	findCandidateAlignments(gw, imageAlignments);
	findCycles(gw, imageAlignments);
}

void BundlerManager::findCycles(GraphWrapper &gw, vector<vector<vector<imageTransformation> > > &imageAlignments)
{
	vector<Cycle> cycles;
	gw.findAll3Cycles(cycles);
	
	string cyclesFilename = mainDirectory + "/matches/cycles.txt";
	ofstream fout_cycles(cyclesFilename.c_str(), ios::out);
	
	EpipolarGeometry eg;
	for(int i=0; i<cycles.size(); i++)
	{
		Matrix3f acc;
		acc.setIdentity();
		for(int m=0; m<cycles[i].path.size(); m++)
		{
			printf("%d ", cycles[i].path[m]);
			fout_cycles << cycles[i].path[m] << " ";
			int img1, img2;
			if(m==cycles[i].path.size()-1)
			{
				img1 = cycles[i].path[m];
				img2 = cycles[i].path[0];
			}
			else
			{
				img1 = cycles[i].path[m];
				img2 = cycles[i].path[m+1];
			}
			if(img1 < img2)
				acc *= imageAlignments[img1][img2][0].fundMatrix;
			else
				acc *= (imageAlignments[img2][img1][0].fundMatrix).transpose();
		}
		printf("\n");
		eg.isIdentityRotation(acc);
		printf("acc:\n");
		printf("%f %f %f\n", acc[0][0], acc[1][0], acc[2][0]);
		printf("%f %f %f\n", acc[0][1], acc[1][1], acc[2][1]);
		printf("%f %f %f\n", acc[0][2], acc[1][2], acc[2][2]);
		fout_cycles << endl;
	}
	
	fout_cycles.close();
}

void BundlerManager::convertAlignmentsToMatches()
{
	EpipolarGeometry epGeo;
	RepetitionFinder rf;
	
	ReadAllFeaturePoints();
	readMatches();
	
	int numImages = imageFileNames.size();
	int noTemplates = repetitions.size();
	
	string alignmentFileName = mainDirectory + "/matches/optGraph.txt";
	ifstream fin(alignmentFileName.c_str(), ios::in);
	
	vector<vector<imageTransformation> > alignments;
	alignments.resize(numImages);
	
	for(int i=0; i<numImages; i++)
		alignments[i].resize(numImages);
	
	GraphWrapper gw(numImages);
	
	if(!fin)
		return ;
	int noAlignments;
	fin >> noAlignments;
	
	vector<int> sourceVertices;
	vector<int> targetVertices;
	
	for(int a=0; a<noAlignments; a++)
	{
		int img1, img2;
		float var, weight;
		int noGridAlg;
		Matrix3f fundMatrix;
		fin >> img1 >> img2 >> var >> weight;
		fin >> noGridAlg;
		printf("*****%d %d %f******\n", img1, img2, var);
		
		imageTransformation tr;
		tr.ind1 = img1; tr.ind2 = img2; tr.score = var;
		
		if(noGridAlg == 0)
			tr.gridMatch = 0;
		else
			tr.gridMatch = 1;
		
		for(int i=0; i<noGridAlg; i++)
		{
			gridShift gs;
			int planeId1, planeId2;
			int gridId1, gridId2;
			int col, row;
			int templateId;
			fin >> templateId >> planeId1 >> planeId2 >> gridId1 >> gridId2 >> col >> row;
			gs.templateId = templateId; gs.planeIndex1 = planeId1; gs.planeIndex2 = planeId2; 
			gs.gridIndex1 = gridId1; gs.gridIndex2 = gridId2;
			gs.shiftCol = col; gs.shiftRow = row;
			tr.gridShifts.push_back(gs);
		}
		fin >> tr.fundMatrix[0][0] >> tr.fundMatrix[1][0] >> tr.fundMatrix[2][0];
		fin >> tr.fundMatrix[0][1] >> tr.fundMatrix[1][1] >> tr.fundMatrix[2][1];
		fin >> tr.fundMatrix[0][2] >> tr.fundMatrix[1][2] >> tr.fundMatrix[2][2];
		if(var > 0.2)
			continue;
		
		gw.addEdge(img1, img2, 1.0/weight);
		
		alignments[img1][img2] = tr;
	}
	
	
	gw.findMinSpanningTree(sourceVertices, targetVertices);
	
	printf("spanning tree:\n");
	for(int i=0;i<sourceVertices.size(); i++)
	{
		printf("%d, %d\n", sourceVertices[i], targetVertices[i]);
		
		if(alignments[sourceVertices[i]][targetVertices[i]].gridMatch)
			continue;
		
		vector<Vec2f> aPts;
		vector<Vec2f> bPts;
		
		for(int m=0; m<matches[sourceVertices[i]][targetVertices[i]-sourceVertices[i]-1].size(); m++)
		{
			Vec2f pos(keyInfo[sourceVertices[i]][matches[sourceVertices[i]][targetVertices[i]-sourceVertices[i]-1][m].m_idx1].x, keyInfo[sourceVertices[i]][matches[sourceVertices[i]][targetVertices[i]-sourceVertices[i]-1][m].m_idx1].y);
			aPts.push_back(pos);
			
			pos = Vec2f(keyInfo[targetVertices[i]][matches[sourceVertices[i]][targetVertices[i]-sourceVertices[i]-1][m].m_idx2].x, keyInfo[targetVertices[i]][matches[sourceVertices[i]][targetVertices[i]-sourceVertices[i]-1][m].m_idx2].y);
			bPts.push_back(pos);
		}
		
		Matrix3f R;
		Vec3f t;
		vector<int> tmpInliers;
		
		int noInliers = epGeo.estimatePose(aPts.size(), aPts, bPts, intrinsicMatrices[sourceVertices[i]], intrinsicMatrices[targetVertices[i]], R, t, tmpInliers);
		
		float ratio = (float) noInliers / (float) aPts.size();
		printf("rotation:\n");
		printf("%f %f %f\n", R[0][0], R[1][0], R[2][0]);
		printf("%f %f %f\n", R[0][1], R[1][1], R[2][1]);
		printf("%f %f %f\n", R[0][2], R[1][2], R[2][2]);
		
		
		bool matchesToGrid = false;
		for(int t=0; t<noTemplates; t++)
		{
			for(int p1=0; p1<noRectifiedImages[sourceVertices[i]]; p1++)
			{
				for(int g1=0; g1<repeatingGrids[0][sourceVertices[i]][p1].size(); g1++)
				{
					for(int p2=0; p2<noRectifiedImages[targetVertices[i]]; p2++)
					{
						for(int g2=0; g2<repeatingGrids[0][targetVertices[i]][p2].size(); g2++)
						{
							Matrix3f Mi = rectifyingHomographies[sourceVertices[i]][p1];
							Mi = Mi.getInverseMatrix();
							Matrix3f Mj = rectifyingHomographies[targetVertices[i]][p2];
							Mj = Mj.getInverseMatrix();
							
							vector<Vec2f> taPts;
							vector<Vec2f> tbPts;
							vector<int> inliers;
							
							for(int m=0; m<tmpInliers.size(); m++)
							{
								inliers.push_back(m);
								
								Vec3f pos = expand2To3(aPts[tmpInliers[m]]);
								pos = Mi*pos;
								pos[0] = pos[0]/pos[2];
								pos[1] = pos[1]/pos[2];
								taPts.push_back(shrink3To2(pos));
								
								pos = expand2To3(bPts[tmpInliers[m]]);
								pos = Mj*pos;
								pos[0] = pos[0]/pos[2];
								pos[1] = pos[1]/pos[2];
								tbPts.push_back(shrink3To2(pos));
							}
							
							float bestScale = 1.0;
							Vec2f bestTranslation = Vec2f(0,0);
							epGeo.estimateRectifiedTransformed(taPts.size(), taPts, tbPts,  bestScale, bestTranslation);
							
							vector<float> bestScales;
							vector<Vec2f> bestTranslations;
							vector<vector<int> > bestInliers;
							bestScales.push_back(bestScale);
							bestTranslations.push_back(bestTranslation);
							bestInliers.push_back(inliers);
							vector<int> shiftCol;
							vector<int> shiftRow;
							vector<vector<int> > tmp;
							int bestCandidate;
							
							float threshold = 50.0;
							if(repeatingGrids[0][targetVertices[i]][p2][g2].xTrans[0] != 0.0 && repeatingGrids[0][targetVertices[i]][p2][g2].xTrans[0]/2.0<threshold)
								threshold = repeatingGrids[0][targetVertices[i]][p2][g2].xTrans[0]/2.0;
							if(repeatingGrids[0][targetVertices[i]][p2][g2].yTrans[1] != 0.0 && repeatingGrids[0][targetVertices[i]][p2][g2].yTrans[1]/2.0<threshold)
								threshold = repeatingGrids[0][targetVertices[i]][p2][g2].yTrans[1]/2.0;
							
							float ratio = rf.translateTransformationToGridRANSAC(repeatingGrids[0][targetVertices[i]][p2][g2], repeatingGrids[0][sourceVertices[i]][p1][g1], taPts, tbPts, bestScales, bestTranslations, bestInliers,
																				 threshold, imageSizes[targetVertices[i]][0], imageSizes[targetVertices[i]][1], shiftCol, shiftRow, tmp, bestCandidate);
							if(bestCandidate > -1 )
							{
								if((float)(tmp[bestCandidate].size()) / (float)(noInliers) > 0.5)
									matchesToGrid = true;
							}
						}
					}
				}
			}
		}
		if(matchesToGrid)
		{
			printf("matches to grid!\n");
		}
	}
	
	GraphWrapper finalSolution(numImages);
	for(int i=0; i<sourceVertices.size(); i++)
	{
		printf("%d %d\n", sourceVertices[i], targetVertices[i]);
		finalSolution.addEdge(sourceVertices[i], targetVertices[i], 1.0);
	}
	
	int maxIter = 3;
	for(int iter = 0; iter<maxIter; iter++ )
	{
		int prevEdges = sourceVertices.size();
		for(int i=0; i<numImages; i++)
		{
			for(int j=i+1; j<numImages; j++)
			{
				edge_descriptor ed;
				if(finalSolution.doesEdgeExist(i, j, ed))
					continue;
			
				for(int k=0; k<numImages; k++)
				{
					edge_descriptor ed1, ed2;
					if(finalSolution.doesEdgeExist(i, k, ed1) && finalSolution.doesEdgeExist(k, j, ed2))
					{
						int col = 0;
						int row = 0;
						int firstPlane, firstGrid;
						int plane, grid;
						int templateId = 0;
					
						bool valid = true;
						if(i < k)
						{
							if(alignments[i][k].gridMatch)
							{
								plane = alignments[i][k].gridShifts[0].planeIndex1;
								grid = alignments[i][k].gridShifts[0].gridIndex1;
							
								firstPlane = alignments[i][k].gridShifts[0].planeIndex1;
								firstGrid = alignments[i][k].gridShifts[0].gridIndex1;
							
								if(noRectifiedImages[i] <= plane)
									valid = false;
								else if(repetitions[templateId][i][plane].size() <= grid)
								{
									valid = false;
								}
								else
								{
									col += alignments[i][k].gridShifts[0].shiftCol;
									row += alignments[i][k].gridShifts[0].shiftRow;
									plane = alignments[i][k].gridShifts[0].planeIndex2;
									grid = alignments[i][k].gridShifts[0].gridIndex2; 
								}
							}
							else
								valid = false;
						}
						else
						{
							if(alignments[k][i].gridMatch)
							{
								plane = alignments[k][i].gridShifts[0].planeIndex2;
								grid = alignments[k][i].gridShifts[0].gridIndex2;
							
								firstPlane = alignments[k][i].gridShifts[0].planeIndex2;
								firstGrid = alignments[k][i].gridShifts[0].gridIndex2;
							
								if(noRectifiedImages[i] <= plane)
									valid = false;
								else if(repetitions[templateId][i][plane].size() <= grid)
								{
									valid = false;
								}
								else
								{
									col -= alignments[k][i].gridShifts[0].shiftCol;
									row -= alignments[k][i].gridShifts[0].shiftRow;
									plane = alignments[k][i].gridShifts[0].planeIndex1;
									grid = alignments[k][i].gridShifts[0].gridIndex1;
								}
							}
							else
								valid = false;
						}
					
						if(!valid)
							continue;
					
						if(k < j)
						{
							if(alignments[k][j].gridMatch)
							{
								if(plane != alignments[k][j].gridShifts[0].planeIndex1 || grid != alignments[k][j].gridShifts[0].gridIndex1)
								{
									valid = false;
								}
								else
								{
									col += alignments[k][j].gridShifts[0].shiftCol;
									row += alignments[k][j].gridShifts[0].shiftRow;
									plane = alignments[k][j].gridShifts[0].planeIndex2;
									grid = alignments[k][j].gridShifts[0].gridIndex2; 
								}
							}
							else
								valid = false;
						}
						else
						{
							if(alignments[j][k].gridMatch)
							{
								if(plane != alignments[j][k].gridShifts[0].planeIndex2 || grid != alignments[j][k].gridShifts[0].gridIndex2)
								{
									valid = false;
								}
								else
								{
									col -= alignments[j][k].gridShifts[0].shiftCol;
									row -= alignments[j][k].gridShifts[0].shiftRow;
									plane = alignments[j][k].gridShifts[0].planeIndex1;
									grid = alignments[j][k].gridShifts[0].gridIndex1;
								}
							}
							else
								valid = false;
						}
					
						if(valid)
						{
							if(noRectifiedImages[j] <= plane)
								valid = false;
							else if(repetitions[templateId][j][plane].size() <= grid)
							{
								valid = false;
							}
							else
							{
								imageTransformation tr;
								tr.gridMatch = true;
								tr.ind1 = i;
								tr.ind2 = j;
								gridShift gs;
								gs.templateId = templateId;
								gs.planeIndex1 = firstPlane; gs.planeIndex2 = plane;
								gs.gridIndex1 = firstGrid; gs.gridIndex2 = grid;
								gs.shiftCol = col; gs.shiftRow = row;
								tr.gridShifts.push_back(gs);
								alignments[i][j] = tr;
							
								printf("new alignment from img %d to img %d\n", i , j);
							
								printf("plane:%d, grid:%d to plane:%d grid:%d:%d, %d\n", firstPlane, firstGrid, plane, grid, col, row);
								sourceVertices.push_back(i); targetVertices.push_back(j);
							
								break;
							}
						}
					}
				
				}
			}
		}
		for(int n=prevEdges; n<sourceVertices.size(); n++)
		{
			finalSolution.addEdge(sourceVertices[n], targetVertices[n], 1.0);
		}
	}
	
	vector<vector<vector<KeypointMatch> > > newMatches;
	newMatches.resize(numImages);
	for(int i=0; i<numImages; i++)
		newMatches[i].resize(numImages-i-1);
	
	vector<registeredImages> groups;
	
	for(int i=0; i<sourceVertices.size(); i++)
	{
		for(int m=0; m<matches[sourceVertices[i]][targetVertices[i]-sourceVertices[i]-1].size(); m++)
		{
			newMatches[sourceVertices[i]][targetVertices[i]-sourceVertices[i]-1].push_back(matches[sourceVertices[i]][targetVertices[i]-sourceVertices[i]-1][m]);
		}
		
		vector<KeypointMatch> tmpMatches;
		imageTransformation tr = alignments[sourceVertices[i]][targetVertices[i]];
		int img1 = tr.ind1;
		int img2 = tr.ind2;
		
		if(img1 == 19 && img2 == 20)
			int debug = 1;
		
		if(tr.gridMatch)
		{
			evaluateImgToImgGridAlignment(img1, img2, tr.gridShifts[0].templateId, tr.gridShifts[0].planeIndex1, tr.gridShifts[0].planeIndex2, 
										  tr.gridShifts[0].gridIndex1, tr.gridShifts[0].gridIndex2, tr.gridShifts[0].shiftCol, tr.gridShifts[0].shiftRow, tmpMatches);
			printf("from img %d to img %d: %d matches\n", img1, img2, tmpMatches.size());
			newMatches[img1][img2-img1-1].clear();
			for(int m=0; m<tmpMatches.size(); m++)
			{
				newMatches[img1][img2-img1-1].push_back(tmpMatches[m]);
			}
			/*for(int m=0; m<matches[img1][img2-img1-1].size(); m++)
				newMatches[img1][img2-img1-1].push_back(matches[img1][img2-img1-1][m]);
             */
		}
		else
		{
			vector<Vec2f> aPts;
			vector<Vec2f> bPts;
			
			for(int m=0; m<matches[img1][img2-img1-1].size(); m++)
			{
				Vec2f pos(keyInfo[img1][matches[img1][img2-img1-1][m].m_idx1].x, keyInfo[img1][matches[img1][img2-img1-1][m].m_idx1].y);
				aPts.push_back(pos);
				
				pos = Vec2f(keyInfo[img2][matches[img1][img2-img1-1][m].m_idx2].x, keyInfo[img2][matches[img1][img2-img1-1][m].m_idx2].y);
				bPts.push_back(pos);
			}
			
			Matrix3f R;
			Vec3f t;
			vector<int> tmpInliers;
			
			int noInliers = epGeo.estimatePose(aPts.size(), aPts, bPts, intrinsicMatrices[sourceVertices[i]], intrinsicMatrices[targetVertices[i]], R, t, tmpInliers);
			newMatches[img1][img2-img1-1].clear();
			//for(int m=0; m<tmpInliers.size(); m++)
			//{
			//	newMatches[img1][img2-img1-1].push_back(matches[img1][img2-img1-1][tmpInliers[m]]);
			//}
			for(int m=0; m<matches[img1][img2-img1-1].size(); m++)
				newMatches[img1][img2-img1-1].push_back(matches[img1][img2-img1-1][m]);
             
		}
		for(int c=0; c<tr.gridShifts.size(); c++)
		{
			bool found = false;
			gridShift gs = tr.gridShifts[c];
			
			for(int g=0; g<groups.size(); g++)
			{
				if(groups[g].templateId != gs.templateId)
					continue;
				vector<int>::iterator it = find(groups[g].imageIndices.begin(), groups[g].imageIndices.end(), tr.ind1);
				if(it == groups[g].imageIndices.end())
					continue;
				int ind = it - groups[g].imageIndices.begin();
				if(gs.planeIndex1 != groups[g].planeIndices[ind])
					continue;
				if(gs.gridIndex1 != groups[g].gridIndices[ind])
					continue;
				//found
				found = true;
				{
					bool mergeGroups = false;
					for(int g2=0; g2<groups.size(); g2++)
					{
						if(groups[g2].templateId != gs.templateId)
							continue;
						vector<int>::iterator it = find(groups[g2].imageIndices.begin(), groups[g2].imageIndices.end(), tr.ind2);
						if(it == groups[g2].imageIndices.end())
							continue;
						int ind2 = it - groups[g2].imageIndices.begin();
						if(gs.planeIndex2 != groups[g2].planeIndices[ind2])
							continue;
						if(gs.gridIndex2 != groups[g2].gridIndices[ind2])
							continue;
						
						mergeGroups = true;
						if(g == g2)
							break;
						
						//merge two groups
						groups[g].imageIndices.push_back(tr.ind2);
						groups[g].planeIndices.push_back(gs.planeIndex2);
						groups[g].gridIndices.push_back(gs.gridIndex2);
						groups[g].shifts.push_back(Vec2i(groups[g].shifts[ind][0]-gs.shiftCol, groups[g].shifts[ind][1]-gs.shiftRow));
						
						int colDiff = groups[g].shifts[groups[g].shifts.size()-1][0] - groups[g2].shifts[ind2][0];
						int rowDiff = groups[g].shifts[groups[g].shifts.size()-1][1] - groups[g2].shifts[ind2][1];
						
						for(int a=0; a<groups[g2].imageIndices.size(); a++)
						{
							if(groups[g2].imageIndices[a] == tr.ind2)
								continue;
							if(find(groups[g].imageIndices.begin(), groups[g].imageIndices.end(), groups[g2].imageIndices[a]) != groups[g].imageIndices.end())
								continue;
							groups[g].imageIndices.push_back(groups[g2].imageIndices[a]);
							groups[g].planeIndices.push_back(groups[g2].planeIndices[a]);
							groups[g].gridIndices.push_back(groups[g2].gridIndices[a]);
							groups[g].shifts.push_back(Vec2i(colDiff+groups[g2].shifts[a][0], rowDiff+groups[g2].shifts[a][1]));
						}
						groups.erase(groups.begin()+g2);
						break;
					}
					if(!mergeGroups)
					{
						groups[g].imageIndices.push_back(tr.ind2);
						groups[g].planeIndices.push_back(gs.planeIndex2);
						groups[g].gridIndices.push_back(gs.gridIndex2);
						groups[g].shifts.push_back(Vec2i(groups[g].shifts[ind][0]-gs.shiftCol, groups[g].shifts[ind][1]-gs.shiftRow));
					}
					break;
				}
			}
			if(!found)
			{
				for(int g=0; g<groups.size(); g++)
				{
					if(groups[g].templateId != gs.templateId)
						continue;
					vector<int>::iterator it = find(groups[g].imageIndices.begin(), groups[g].imageIndices.end(), tr.ind2);
					if(it == groups[g].imageIndices.end())
						continue;
					int ind = it - groups[g].imageIndices.begin();
					if(gs.planeIndex2 != groups[g].planeIndices[ind])
						continue;
					if(gs.gridIndex2 != groups[g].gridIndices[ind])
						continue;
					//found
					found = true;
					{
						bool mergeGroups = false;
						for(int g2=0; g2<groups.size(); g2++)
						{
							if(groups[g2].templateId != gs.templateId)
								continue;
							vector<int>::iterator it = find(groups[g2].imageIndices.begin(), groups[g2].imageIndices.end(), tr.ind1);
							if(it == groups[g2].imageIndices.end())
								continue;
							int ind2 = it - groups[g2].imageIndices.begin();
							if(gs.planeIndex1 != groups[g2].planeIndices[ind])
								continue;
							if(gs.gridIndex1 != groups[g2].gridIndices[ind])
								continue;
						
							mergeGroups = true;
							if(g == g2)
								break;
						
							//merge two groups
							groups[g].imageIndices.push_back(tr.ind1);
							groups[g].planeIndices.push_back(gs.planeIndex1);
							groups[g].gridIndices.push_back(gs.gridIndex1);
							groups[g].shifts.push_back(Vec2i(groups[g].shifts[ind][0]+gs.shiftCol, groups[g].shifts[ind][1]+gs.shiftRow));
						
							int colDiff = groups[g].shifts[groups[g].shifts.size()-1][0] - groups[g2].shifts[ind2][0];
							int rowDiff = groups[g].shifts[groups[g].shifts.size()-1][1] - groups[g2].shifts[ind2][1];
						
							for(int a=0; a<groups[g2].imageIndices.size(); a++)
							{
								if(groups[g2].imageIndices[a] == tr.ind1)
									continue;
								if(find(groups[g].imageIndices.begin(), groups[g].imageIndices.end(), groups[g2].imageIndices[a]) != groups[g].imageIndices.end())
									continue;
								groups[g].imageIndices.push_back(groups[g2].imageIndices[a]);
								groups[g].planeIndices.push_back(groups[g2].planeIndices[a]);
								groups[g].gridIndices.push_back(groups[g2].gridIndices[a]);
								groups[g].shifts.push_back(Vec2i(colDiff+groups[g2].shifts[a][0], rowDiff+groups[g2].shifts[a][1]));
							}
							groups.erase(groups.begin()+g2);
							break;
						}
						if(!mergeGroups)
						{
							groups[g].imageIndices.push_back(tr.ind1);
							groups[g].planeIndices.push_back(gs.planeIndex1);
							groups[g].gridIndices.push_back(gs.gridIndex1);
							groups[g].shifts.push_back(Vec2i(groups[g].shifts[ind][0]+gs.shiftCol, groups[g].shifts[ind][1]+gs.shiftRow));
						}
						break;
					}
				}
			}
			if(!found)
			{
				registeredImages ri;
				ri.templateId = gs.templateId;
				ri.imageIndices.push_back(tr.ind1); ri.imageIndices.push_back(tr.ind2);
				ri.planeIndices.push_back(gs.planeIndex1); ri.planeIndices.push_back(gs.planeIndex2);
				ri.gridIndices.push_back(gs.gridIndex1); ri.gridIndices.push_back(gs.gridIndex2);
				ri.shifts.push_back(Vec2i(0,0));
				ri.shifts.push_back(Vec2i(-gs.shiftCol, -gs.shiftRow));
				groups.push_back(ri);
			}
		}
	}
	
	//int iter = find(groups[1].imageIndices.begin(), groups[1].imageIndices.end(), 22) - groups[1].imageIndices.begin();
	//groups[1].imageIndices.push_back(19); groups[1].planeIndices.push_back(0); groups[1].gridIndices.push_back(0);
	//groups[1].shifts.push_back(groups[1].shifts[iter]-Vec2i(0, 0));
	
	//groups[1].imageIndices.push_back(20); groups[1].planeIndices.push_back(0); groups[1].gridIndices.push_back(0);
	//groups[1].shifts.push_back(groups[1].shifts[iter]-Vec2i(0, 0));
	
	//groups.erase(groups.begin()+2);
	
	vector<repetitionGrid> finalGrids;
	vector<vector<int> > groupIndices;
	vector<int> gridSize;
	vector<vector<vector<vector<Vec2f> > > > gridCoord;
	vector<vector<vector<vector<int> > > > gridCoordPlanes;
	vector<vector<vector<vector<int> > > > gridCoordGrids;
	
	int totalNoPts = 0;
	
	for(int g=0; g<groups.size(); g++)
	{
		vector<int> curGroupIndices;
		vector<vector<vector<Vec2f> > > curGridCoord;
		vector<vector<vector<int> > > curGridCoordPlanes;
		vector<vector<vector<int> > > curGridCoordGrids;
		
		curGridCoord.resize(numImages);
		curGridCoordPlanes.resize(numImages);
		curGridCoordGrids.resize(numImages);
		curGroupIndices.resize(numImages, -1);
		
		int minCol, minRow;
		for (int i=0; i<groups[g].imageIndices.size(); i++) 
		{
			curGroupIndices[groups[g].imageIndices[i]] = g;
			if(i==0)
			{
				minCol = groups[g].shifts[i][0];
				minRow = groups[g].shifts[i][1];
			}
			else
			{
				if(groups[g].shifts[i][0] < minCol)
					minCol = groups[g].shifts[i][0];
				if(groups[g].shifts[i][1] < minRow)
					minRow = groups[g].shifts[i][1];
			}
		}
		
		printf("******* group ******\n");
		int startIndex = -1;
		int selCol, selRow;
		for (int i=0; i<groups[g].imageIndices.size(); i++) 
		{
			groups[g].shifts[i][0] = groups[g].shifts[i][0] - minCol;
			groups[g].shifts[i][1] = groups[g].shifts[i][1] - minRow;
			
			if(startIndex == -1)
			{
				selCol = groups[g].shifts[i][0];
				selRow = groups[g].shifts[i][1];
				startIndex = i;
			}
			else if(groups[g].shifts[i][0] < selCol || groups[g].shifts[i][1] < selRow)
			{
				selCol = groups[g].shifts[i][0];
				selRow = groups[g].shifts[i][1];
				startIndex = i;
			}
			printf("img %d, plane %d, grid %d:%d, %d\n", groups[g].imageIndices[i], groups[g].planeIndices[i], groups[g].gridIndices[i], groups[g].shifts[i][0], groups[g].shifts[i][1]);
		}
		
		if(startIndex == -1)
		{
			printf("Error: Grid merging not possible!\n");
			continue;
		}
		repetitionGrid finalGrid = repeatingGrids[groups[g].templateId][groups[g].imageIndices[startIndex]][groups[g].planeIndices[startIndex]][groups[g].gridIndices[startIndex]];
		for(int i=0; i<groups[g].imageIndices.size(); i++)
		{
			if(i==startIndex)
				continue;
			finalGrid = rf.mergeGrids(finalGrid, repeatingGrids[groups[g].templateId][groups[g].imageIndices[i]][groups[g].planeIndices[i]][groups[g].gridIndices[i]], 
									  groups[g].shifts[i][0], groups[g].shifts[i][1]);
		}
		vector<vector<vector<int> > > visImageIndices;
		int noCols = finalGrid.numColumns;
		int noRows = finalGrid.numRows;
		visImageIndices.resize(noCols);
		for(int c=0; c<noCols; c++)
			visImageIndices[c].resize(noRows);
	
		for(int i=0; i<groups[g].imageIndices.size(); i++)
		{
			curGridCoord[groups[g].imageIndices[i]].resize(noCols);
			curGridCoordPlanes[groups[g].imageIndices[i]].resize(noCols);
			curGridCoordGrids[groups[g].imageIndices[i]].resize(noCols);
			for(int c=0; c<noCols; c++)
			{
				curGridCoord[groups[g].imageIndices[i]][c].resize(noRows, Vec2f(0.0, 0.0));
				curGridCoordPlanes[groups[g].imageIndices[i]][c].resize(noRows, -1);
				curGridCoordGrids[groups[g].imageIndices[i]][c].resize(noRows, -1);
			}
			
			repetitionGrid grid = repeatingGrids[groups[g].templateId][groups[g].imageIndices[i]][groups[g].planeIndices[i]][groups[g].gridIndices[i]];
			for(int c=0; c<grid.gridCells.size(); c++)
			{
				visImageIndices[grid.gridCells[c].colNumber+groups[g].shifts[i][0]][grid.gridCells[c].rowNumber+groups[g].shifts[i][1]].push_back(groups[g].imageIndices[i]);
				if(grid.gridCells[c].temporary)
					printf("temporary\n");
				if(!grid.gridCells[c].temporary)
				{
					curGridCoord[groups[g].imageIndices[i]][grid.gridCells[c].colNumber+groups[g].shifts[i][0]][grid.gridCells[c].rowNumber+groups[g].shifts[i][1]] = grid.gridCells[c].center;
					curGridCoordPlanes[groups[g].imageIndices[i]][grid.gridCells[c].colNumber+groups[g].shifts[i][0]][grid.gridCells[c].rowNumber+groups[g].shifts[i][1]] = groups[g].planeIndices[i];
					curGridCoordGrids[groups[g].imageIndices[i]][grid.gridCells[c].colNumber+groups[g].shifts[i][0]][grid.gridCells[c].rowNumber+groups[g].shifts[i][1]] = groups[g].gridIndices[i];
				}
			}
		}
		printf("final grid\n");
		for(int c=0; c<noCols; c++)
		{
			for(int r=0; r<noRows; r++)
			{
				printf("cell at col %d, row %d is visible in:", c, r);
				for (int i=0; i<visImageIndices[c][r].size(); i++) 
				{
					printf("%d ", visImageIndices[c][r][i]);
				}
				printf("\n");
			}
		}
		
		finalGrids.push_back(finalGrid);
		gridSize.push_back(finalGrid.numColumns * finalGrid.numRows);
		totalNoPts += (finalGrid.numColumns * finalGrid.numRows);
		
		for(int i=0; i<numImages; i++)
		{
			if(curGroupIndices[i] == -1)
				curGroupIndices[i] = finalGrids.size();
		}
		
		groupIndices.push_back(curGroupIndices);
		gridCoord.push_back(curGridCoord);
		gridCoordPlanes.push_back(curGridCoordPlanes);
		gridCoordGrids.push_back(curGridCoordGrids);
	}
	
	string gridFileName = mainDirectory + "/matches/repetitions.txt";
	ofstream gridFout(gridFileName.c_str(), ios::out);
	string gridCellsFileName = mainDirectory + "/matches/repetitionCells.txt";
	ofstream gridCellsFout(gridCellsFileName.c_str(), ios::out);
	gridFout << numImages << " " << finalGrids.size() << endl;
	gridCellsFout << numImages << " " << finalGrids.size() << endl;
	
	for(int g=0; g<finalGrids.size(); g++)
	{
		gridFout << finalGrids[g].numRows << " " << finalGrids[g].numColumns << endl;
		gridCellsFout << finalGrids[g].numRows << " " << finalGrids[g].numColumns << endl;
		
		for(int i=0; i<numImages; i++)
		{
			printf("image %d\n", i);
			if(groupIndices[g][i] == g)
			{
				for(int r=0; r<finalGrids[g].numRows; r++)
				{
					for(int c=0; c<finalGrids[g].numColumns; c++)
					{
						Vec2f pos = gridCoord[g][i][c][r];
						gridCellsFout << gridCoordPlanes[g][i][c][r] << " " << gridCoordGrids[g][i][c][r] << " " << pos[0] << " " << pos[1] << " ";
						printf("r:%d c:%d: %f %f\n", r, c, pos[0], pos[1]);
						if(pos == Vec2f(0.0, 0.0))
							gridFout << "-1 ";
						else
							gridFout << "1 ";
					}
				}
			}
			else
			{
				for(int r=0; r<finalGrids[g].numRows; r++)
				{
					for(int c=0; c<finalGrids[g].numColumns; c++)
					{
						gridFout << "-1 ";
						gridCellsFout << "-1 -1 0.0 0.0 ";
					}
				}
			}
			gridFout << "\n";
			gridCellsFout << "\n";
		}
	}
	gridFout.close();
	gridCellsFout.close();
	
	for(int i=0; i<numImages; i++)
	{
		char keyFile[1024];
		sprintf(keyFile, "%08d.%s", i, "key");
		string filename = mainDirectory + "/matches/" + keyFile;
		ofstream fout(filename.c_str());
		
		fout << keyInfo[i].size() + totalNoPts << " 128" << endl;
		
		for(int g=0; g<finalGrids.size(); g++)
		{
			if(groupIndices[g][i] == g)
			{
				int planeIndex = groups[g].planeIndices[find(groups[g].imageIndices.begin(), groups[g].imageIndices.end(), i) - groups[g].imageIndices.begin()];
				for(int r=0; r<finalGrids[g].numRows; r++)
				{
					for(int c=0; c<finalGrids[g].numColumns; c++)
					{
						Vec3f pos = Vec3f(gridCoord[g][i][c][r][0], gridCoord[g][i][c][r][1], 1.0);
						if(pos != Vec3f(0.0, 0.0, 1.0))
						{
							pos = rectifyingHomographies[i][planeIndex] * pos;
							pos[0] /= pos[2];
							pos[1] /= pos[2];
						}
						
						fout << pos[1] << " " << pos[0] << " 0.0 0.0 1\n";
						for(int d=0; d<7; d++)
						{
							if(d==6)
								fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
							else
								fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
						}
					}
				}
			}
			else
			{
				for(int c=0; c<gridSize[g]; c++)
				{
					fout << "0.0 0.0 0.0 0.0 1\n";
					for(int d=0; d<7; d++)
					{
						if(d==6)
							fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
						else
							fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
					}
				}
			}
		}
		
		/*for(int g=0; g<groupIndices[g][i]; g++)
		{
			for(int c=0; c<gridSize[g]; c++)
			{
				fout << "0.0 0.0 0.0 0.0 1\n";
				for(int d=0; d<7; d++)
				{
					if(d==6)
						fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
					else
						fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
				}
			}
		}
		if(groupIndices[g][i] < finalGrids.size())
		{
			int planeIndex = groups[groupIndices[g][i]].planeIndices[find(groups[groupIndices[g][i]].imageIndices.begin(), groups[groupIndices[g][i]].imageIndices.end(), i) - groups[groupIndices[g][i]].imageIndices.begin()];
			for(int r=0; r<finalGrids[groupIndices[g][i]].numRows; r++)
			{
				for(int c=0; c<finalGrids[groupIndices[g][i]].numColumns; c++)
				{
					
					Vec3f pos = rectifyingHomographies[i][planeIndex] * Vec3f(gridCoord[g][i][c][r][0], gridCoord[g][i][c][r][1], 1.0);
					pos[0] /= pos[2];
					pos[1] /= pos[2];
					
					fout << pos[1] << " " << pos[0] << " 0.0 0.0 1\n";
					for(int d=0; d<7; d++)
					{
						if(d==6)
							fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
						else
							fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
					}
				}
			}
		}
		for(int g=groupIndices[g][i]+1; g<finalGrids.size(); g++)
		{
			for(int c=0; c<gridSize[g]; c++)
			{
				fout << "0.0 0.0 0.0 0.0 1\n";
				for(int d=0; d<7; d++)
				{
					if(d==6)
						fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
					else
						fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
				}
			}
		}*/
		
		for(int r=0; r<keyInfo[i].size(); r++)
		{
			fout << keyInfo[i][r].y << " " << keyInfo[i][r].x << " 0.0 0.0 0\n";
			for(int d=0; d<7; d++)
			{
				if(d==6)
					fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
				else
					fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
			}
		}
		
		fout.close();
	}
	
	for(int i=0; i<numImages; i++)
	{
		for(int j=i+1; j<numImages; j++)
		{
			matches[i][j-i-1].clear();
			for(int g=0; g<finalGrids.size(); g++)
			{
				if(groupIndices[g][i] == groupIndices[g][j] && groupIndices[g][i] == g)
				{
					int noPrevPts = 0;
					for(int gp=0; gp<g; gp++)
						noPrevPts += gridSize[gp];
				
					//grid matches
					for(int r=0; r<finalGrids[g].numRows; r++)
					{
						for(int c=0; c<finalGrids[g].numColumns; c++)
						{
							Vec2f pos1 = gridCoord[g][i][c][r];
							Vec2f pos2 = gridCoord[g][j][c][r];
							if(pos1 != Vec2f(0.0, 0.0) && pos2 != Vec2f(0.0, 0.0))
							{
								KeypointMatch kp;
								kp.m_idx1 = noPrevPts + r*finalGrids[g].numColumns + c;
								kp.m_idx2 = kp.m_idx1; 
								matches[i][j-i-1].push_back(kp);
							}
						}
					}
				}
			}
			
			for(int m=0; m<newMatches[i][j-i-1].size(); m++)
			{
				KeypointMatch kp = newMatches[i][j-i-1][m];
				kp.m_idx1 = kp.m_idx1 + totalNoPts;
				kp.m_idx2 = kp.m_idx2 + totalNoPts;
				matches[i][j-i-1].push_back(kp);
			}
		}
	}
}

void BundlerManager::findCandidateAlignments(GraphWrapper &gw, vector<vector<vector<imageTransformation> > > &imageAlignments)
{
	EpipolarGeometry epGeo;
	RepetitionFinder rf;
	
	int numImages = imageFileNames.size();
	int noTemplates = repetitions.size();
	
	int minNoMatches = 12;
	float graphEdgeThreshold = 0.5;
	
	//vector<vector<vector<imageTransformation> > > imageAlignments;
	imageAlignments.resize(numImages);
	
	string alignmentFileName = mainDirectory + "/matches/graph.txt";
	int alignmentCount = 0;
	int maxScore = 0;
	
	ReadAllFeaturePoints();
	readMatches();
	
	ifstream fin(alignmentFileName.c_str());
	if(fin)
	{
		for(int i=0; i<numImages; i++)
		{
			imageAlignments[i].resize(numImages);
		}
		
		fin >> alignmentCount;
		for(int a =0; a<alignmentCount; a++)
		{
			imageTransformation tr;
			fin >> tr.ind1 >> tr.ind2 >> tr.score >> tr.gridMatch;
			fin >> tr.scale >> tr.translation[0] >> tr.translation[1];
			int noGridShifts;
			fin >> noGridShifts;
			for(int g=0; g<noGridShifts; g++)
			{
				gridShift gs;
				fin >> gs.templateId >> gs.planeIndex1 >> gs.planeIndex2 >> gs.gridIndex1 >> gs.gridIndex2 >> gs.shiftCol >> gs.shiftRow;
				gs.imageIndex1 = tr.ind1; gs.imageIndex2 = tr.ind2;
				tr.gridShifts.push_back(gs);
			}
			fin >> tr.fundMatrix[0][0] >> tr.fundMatrix[1][0] >> tr.fundMatrix[2][0];
			fin >> tr.fundMatrix[0][1] >> tr.fundMatrix[1][1] >> tr.fundMatrix[2][1];
			fin >> tr.fundMatrix[0][2] >> tr.fundMatrix[1][2] >> tr.fundMatrix[2][2];
			
			imageAlignments[tr.ind1][tr.ind2].push_back(tr);
			gw.addEdge(tr.ind1, tr.ind2, tr.score);
		}
		
		fin.close();
		
		return;
	}
	
	for(int i=0; i<numImages; i++)
	{
		imageAlignments[i].resize(numImages);
		for(int j=0; j<numImages; j++)
			imageAlignments[i][j].resize(1);
	}
	
	for(int i=0; i<numImages; i++)
	{
		
		for(int j=i+1; j<numImages; j++)
		{
			printf("********* img %d to img %d ***********\n", i, j);
			int maxNoMatches = 0;
			Vec2i bestAlignment;
			int templateId;
			int planeIndex1, planeIndex2;
			int grid1, grid2;
			vector<KeypointMatch> bestMatches;
			
			bool gridFound = false;
			
			imageAlignments[i][j][0].ind1 = i; imageAlignments[i][j][0].ind2 = j;
			
			for(int t=0; t<noTemplates; t++)
			{
				for(int p1=0; p1<noRectifiedImages[i]; p1++)
				{
					for(int g1=0; g1<repeatingGrids[t][i][p1].size(); g1++)
					{
						for(int p2=0; p2<noRectifiedImages[j]; p2++)
						{
							for(int g2=0; g2<repeatingGrids[t][j][p2].size(); g2++)
							{
								gridFound = true;
								int c, r, matchSize;
								vector<KeypointMatch> tmpMatches;
								evaluateAllPossibleAlignments(i, j, t, p1, p2, g1, g2, c, r, matchSize, tmpMatches);
								if(matchSize > maxNoMatches)
								{
									maxNoMatches = matchSize;
									templateId = t;
									bestAlignment[0] = c; bestAlignment[1] = r;
									planeIndex1 = p1; planeIndex2 = p2;
									grid1 = g1; grid2 = g2;
									
									bestMatches.clear();
									for(int m=0; m<tmpMatches.size(); m++)
										bestMatches.push_back(tmpMatches[m]);
								}
							}
						}
					}
				}
			}
			
			printf("maxNoMatches:%d\n", maxNoMatches);
			if(maxNoMatches > minNoMatches)
			{
				imageAlignments[i][j][0].gridMatch = true;
				vector<gridShift> grids;
				
				if(maxNoMatches > maxScore)
				{
					maxScore = maxNoMatches;
				}
				
				gridShift gs;
				gs.imageIndex1 = i;
				gs.imageIndex2 = j;
				gs.templateId = templateId;
				gs.planeIndex1 = planeIndex1;
				gs.planeIndex2 = planeIndex2;
				gs.gridIndex1 = grid1;
				gs.gridIndex2 = grid2;
				gs.shiftCol = bestAlignment[0];
				gs.shiftRow = bestAlignment[1];
				imageAlignments[i][j][0].gridShifts.push_back(gs);
				
				printf("best alignment from grid %d on plane %d to grid %d on plane %d\n", grid1, planeIndex1, grid2, planeIndex1);
				
				//change this to alignments between other grids
				Matrix3f Mi = rectifyingHomographies[i][planeIndex1];
				Mi = Mi.getInverseMatrix();
				Matrix3f Mj = rectifyingHomographies[j][planeIndex2];
				Mj = Mj.getInverseMatrix();
				
				vector<Vec2f> aPts;
				vector<Vec2f> bPts;
				vector<int> inliers;
				
				vector<Vec2f> aNrmPts;
				vector<Vec2f> bNrmPts;
				
				for(int m=0; m<bestMatches.size(); m++)
				{
					Vec3f pos(keyInfo[i][bestMatches[m].m_idx1].x, keyInfo[i][bestMatches[m].m_idx1].y, 1.0);
					aNrmPts.push_back(shrink3To2(pos));
					pos = Mi*pos;
					pos[0] /= pos[2];
					pos[1] /= pos[2];
					aPts.push_back(shrink3To2(pos));
					
					pos = Vec3f(keyInfo[j][bestMatches[m].m_idx2].x, keyInfo[j][bestMatches[m].m_idx2].y, 1.0);
					bNrmPts.push_back(shrink3To2(pos));
					pos = Mj*pos;
					pos[0] /= pos[2];
					pos[1] /= pos[2];
					bPts.push_back(shrink3To2(pos));
					
					inliers.push_back(m);
				}
				
				Matrix3f R;
				Vec3f t;
				vector<int> tmpInliers;
				epGeo.estimatePose(aNrmPts.size(), aNrmPts, bNrmPts, intrinsicMatrices[i], intrinsicMatrices[j], R, t, tmpInliers);
				imageAlignments[i][j][0].fundMatrix = R;
				
				float bestScale = 1.0;
				Vec2f bestTranslation = Vec2f(0,0);
				epGeo.estimateRectifiedTransformed(aPts.size(), aPts, bPts,  bestScale, bestTranslation);
				imageAlignments[i][j][0].scale = bestScale;
				imageAlignments[i][j][0].translation = bestTranslation;
				imageAlignments[i][j][0].score = aPts.size();
				
				vector<float> bestScales;
				vector<Vec2f> bestTranslations;
				vector<vector<int> > bestInliers;
				bestScales.push_back(bestScale);
				bestTranslations.push_back(bestTranslation);
				bestInliers.push_back(inliers);
				
				bool valid = true;
				
				for(int t=0; t<noTemplates; t++)
				{
					for(int p1=0; p1<noRectifiedImages[i]; p1++)
					{
						for(int g1=0; g1<repeatingGrids[t][i][p1].size(); g1++)
						{
							if(p1 == planeIndex1 && g1 == grid1 && t == templateId)
								continue;
							for(int p2=0; p2<noRectifiedImages[j]; p2++)
							{
								for(int g2=0; g2<repeatingGrids[t][j][p2].size(); g2++)
								{
									if(p2 == planeIndex2 && g2 == grid2 && t==templateId)
										continue;
									vector<int> shiftCol;
									vector<int> shiftRow;
									vector<vector<int> > tmpInliers;
									int bestCandidate;
									
									float threshold = 50.0;
									if(repeatingGrids[t][j][p2][g2].xTrans[0] != 0.0 && repeatingGrids[t][j][p2][g2].xTrans[0]/2.0<threshold)
										threshold = repeatingGrids[t][j][p2][g2].xTrans[0]/2.0;
									if(repeatingGrids[t][j][p2][g2].yTrans[1] != 0.0 && repeatingGrids[t][j][p2][g2].yTrans[1]/2.0<threshold)
										threshold = repeatingGrids[t][j][p2][g2].yTrans[1]/2.0;
									
									float ratio = rf.translateTransformationToGridRANSAC(repeatingGrids[t][j][p2][g2], repeatingGrids[t][i][p1][g1], aPts, bPts, bestScales, bestTranslations, bestInliers,
																						 threshold, imageSizes[j][0], imageSizes[j][1], shiftCol, shiftRow, tmpInliers, bestCandidate);
									if(shiftCol.size() > 0 && (ratio > graphEdgeThreshold || tmpInliers[bestCandidate].size() > minNoMatches))
									{
										gridShift gs;
										gs.imageIndex1 = i;
										gs.imageIndex2 = j;
										gs.templateId = t;
										gs.planeIndex1 = p1;
										gs.planeIndex2 = p2;
										gs.gridIndex1 = g1;
										gs.gridIndex2 = g2;
										gs.shiftCol = shiftCol[bestCandidate];
										gs.shiftRow = shiftRow[bestCandidate];
										imageAlignments[i][j][0].gridShifts.push_back(gs);
										
										printf("alignment from grid %d on plane %d to grid %d on plane %d\n", g1, p1, g2, p2);
									}
									else
									{
										valid = false;
									}
								}
							}
						}
					}
				}
				if(valid)
				{
					imageAlignments[i][j][0].addedToGraph = true;
					alignmentCount += 1;
				}
				else
				{
					imageAlignments[i][j][0].addedToGraph = false;
					imageAlignments[i][j][0].gridShifts.clear();
					imageAlignments[i][j][0].score = 0.0;
				}
			}
			//else
			{
				/*if(gridFound && noRectifiedImages[i] == 1 && noRectifiedImages[j] == 1)
				{
					it.gridMatch = true;
					it.score = 0.0;
					it.addedToGraph = false;
				}
				else*/
				{
					vector<Vec2f> aPts;
					vector<Vec2f> bPts;
					
					for(int m=0; m<matches[i][j-i-1].size(); m++)
					{
						Vec2f pos(keyInfo[i][matches[i][j-i-1][m].m_idx1].x, keyInfo[i][matches[i][j-i-1][m].m_idx1].y);
						aPts.push_back(pos);
						
						pos = Vec2f(keyInfo[j][matches[i][j-i-1][m].m_idx2].x, keyInfo[j][matches[i][j-i-1][m].m_idx2].y);
						bPts.push_back(pos);
					}
					
					Matrix3f R;
					Vec3f t;
					vector<int> tmpInliers;
					
					int noInliers = epGeo.estimatePose(aPts.size(), aPts, bPts, intrinsicMatrices[i], intrinsicMatrices[j], R, t, tmpInliers);
					
					float ratio = (float) noInliers / (float) aPts.size();
					printf("rotation:\n");
					printf("%f %f %f\n", R[0][0], R[1][0], R[2][0]);
					printf("%f %f %f\n", R[0][1], R[1][1], R[2][1]);
					printf("%f %f %f\n", R[0][2], R[1][2], R[2][2]);
					
					
					bool matchesToGrid = false;
					for(int t=0; t<noTemplates; t++)
					{
						for(int p1=0; p1<noRectifiedImages[i]; p1++)
						{
							for(int g1=0; g1<repeatingGrids[t][i][p1].size(); g1++)
							{
								for(int p2=0; p2<noRectifiedImages[j]; p2++)
								{
									for(int g2=0; g2<repeatingGrids[t][j][p2].size(); g2++)
									{
										Matrix3f Mi = rectifyingHomographies[i][p1];
										Mi = Mi.getInverseMatrix();
										Matrix3f Mj = rectifyingHomographies[j][p2];
										Mj = Mj.getInverseMatrix();
										
										vector<Vec2f> taPts;
										vector<Vec2f> tbPts;
										vector<int> inliers;
										
										for(int m=0; m<tmpInliers.size(); m++)
										{
											inliers.push_back(m);
											
											Vec3f pos = expand2To3(aPts[tmpInliers[m]]);
											pos = Mi*pos;
											pos[0] = pos[0]/pos[2];
											pos[1] = pos[1]/pos[2];
											taPts.push_back(shrink3To2(pos));
											
											pos = expand2To3(bPts[tmpInliers[m]]);
											pos = Mj*pos;
											pos[0] = pos[0]/pos[2];
											pos[1] = pos[1]/pos[2];
											tbPts.push_back(shrink3To2(pos));
										}
										
										float bestScale = 1.0;
										Vec2f bestTranslation = Vec2f(0,0);
										epGeo.estimateRectifiedTransformed(taPts.size(), taPts, tbPts,  bestScale, bestTranslation);
										
										vector<float> bestScales;
										vector<Vec2f> bestTranslations;
										vector<vector<int> > bestInliers;
										bestScales.push_back(bestScale);
										bestTranslations.push_back(bestTranslation);
										bestInliers.push_back(inliers);
										vector<int> shiftCol;
										vector<int> shiftRow;
										vector<vector<int> > tmp;
										int bestCandidate;
										
										float threshold = 50.0;
										if(repeatingGrids[t][j][p2][g2].xTrans[0] != 0.0 && repeatingGrids[t][j][p2][g2].xTrans[0]/2.0<threshold)
											threshold = repeatingGrids[t][j][p2][g2].xTrans[0]/2.0;
										if(repeatingGrids[t][j][p2][g2].yTrans[1] != 0.0 && repeatingGrids[t][j][p2][g2].yTrans[1]/2.0<threshold)
											threshold = repeatingGrids[t][j][p2][g2].yTrans[1]/2.0;
										
										float ratio = rf.translateTransformationToGridRANSAC(repeatingGrids[t][j][p2][g2], repeatingGrids[t][i][p1][g1], taPts, tbPts, bestScales, bestTranslations, bestInliers,
																							 threshold, imageSizes[j][0], imageSizes[j][1], shiftCol, shiftRow, tmp, bestCandidate);
										if(bestCandidate > -1 )
										{
											if((float)(tmp[bestCandidate].size()) / (float)(noInliers) > 0.5)
												matchesToGrid = true;
										}
									}
								}
							}
						}
					}
					
					if(maxNoMatches > minNoMatches && noInliers > maxNoMatches && !matchesToGrid)
					{
						imageAlignments[i][j][0].gridMatch = false;
						imageAlignments[i][j][0].fundMatrix = R;
						imageAlignments[i][j][0].score = noInliers;
						imageAlignments[i][j][0].addedToGraph = true;
						alignmentCount += 1;
						if(noInliers > maxScore)
							maxScore = noInliers;
						printf("noMatches:%d\n", noInliers);
						printf("added to graph without grids\n");
					}
					else if((/*ratio > graphEdgeThreshold ||*/ noInliers > minNoMatches ) && !matchesToGrid )
					{
						imageAlignments[i][j][0].gridMatch = false;
						imageAlignments[i][j][0].fundMatrix = R;
						imageAlignments[i][j][0].score = noInliers;
						imageAlignments[i][j][0].addedToGraph = true;
						alignmentCount += 1;
						if(noInliers > maxScore)
							maxScore = noInliers;
						printf("noMatches:%d\n", noInliers);
						printf("added to graph without grids\n");
					}
					else if(!matchesToGrid)
					{
						imageAlignments[i][j][0].gridMatch = false;
						imageAlignments[i][j][0].fundMatrix = R;
						imageAlignments[i][j][0].score = 0.0;
						imageAlignments[i][j][0].addedToGraph = false;
					}
				}
			}
		}
	}
	
	ofstream alignOut(alignmentFileName.c_str(), ios::out);
	alignOut << alignmentCount << endl;
	for(int i=0; i<numImages; i++)
	{
		for(int j=i+1; j<numImages; j++)
		{
			for(int c=0; c<imageAlignments[i][j].size(); c++)
			{
				//write to file
				imageTransformation tr = imageAlignments[i][j][c];
				tr.score = tr.score / (float) maxScore;
				imageAlignments[i][j][c].score = tr.score;
				if(tr.addedToGraph && tr.ind1 != -1 && tr.ind2 != -1)
				{
					gw.addEdge(i, j, tr.score);
					
					alignOut << tr.ind1 << " " << tr.ind2 << " " << tr.score << " ";
					if(tr.gridMatch)
						alignOut << "1\n";
					else
						alignOut << "0\n";
					
					alignOut << tr.scale << " " << tr.translation[0] << " " << tr.translation[1] << endl;
					alignOut << tr.gridShifts.size() << endl;
					for(int s=0; s<tr.gridShifts.size(); s++)
					{
						alignOut << tr.gridShifts[s].templateId << " " << tr.gridShifts[s].planeIndex1 << " " << tr.gridShifts[s].planeIndex2 << " ";
						alignOut << tr.gridShifts[s].gridIndex1 << " " << tr.gridShifts[s].gridIndex2 << " ";
						alignOut << tr.gridShifts[s].shiftCol << " " << tr.gridShifts[s].shiftRow << endl;
					}
					alignOut << tr.fundMatrix[0][0] << " " << tr.fundMatrix[1][0] << " " << tr.fundMatrix[2][0] << endl;
					alignOut << tr.fundMatrix[0][1] << " " << tr.fundMatrix[1][1] << " " << tr.fundMatrix[2][1] << endl;
					alignOut << tr.fundMatrix[0][2] << " " << tr.fundMatrix[1][2] << " " << tr.fundMatrix[2][2] << endl;
				}
			}
		}
	}
	alignOut.close();
}

void BundlerManager::evaluateAllPossibleAlignments(int refIndex, int neighIndex, int templateId, int refPlaneIndex, int neighPlaneIndex, int refGridIndex, int neighGridIndex, 
												   int &col, int &row, int &matchSize, vector<KeypointMatch> &bestMatches)
{
	RepetitionFinder rf;
	EpipolarGeometry epGeo;
	
	rf.fillMissingGridCellsTemporarily(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex], repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[0].size);
	int minColRef, maxColRef, minRowRef, maxRowRef;
	
	for(int c=0; c<repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells.size(); c++)
	{
		if(c==0)
		{
			minColRef = maxColRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			minRowRef = maxRowRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
		}
		else
		{
			if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber < minColRef)
				minColRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			else if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber > maxColRef)
				maxColRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			
			if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber < minRowRef)
				minRowRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
			else if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber > maxRowRef)
				maxRowRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
		}
	}
	
	rf.fillMissingGridCellsTemporarily(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex], repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[0].size);
	int minColNeigh, maxColNeigh, minRowNeigh, maxRowNeigh;
	
	for(int c=0; c<repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells.size(); c++)
	{
		if(c==0)
		{
			minColNeigh = maxColNeigh = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			minRowNeigh = maxRowNeigh = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
		}
		else
		{
			if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber < minColNeigh)
				minColNeigh = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			else if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber > maxColNeigh)
				maxColNeigh = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			
			if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber < minRowNeigh)
				minRowNeigh = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
			else if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber > maxRowNeigh)
				maxRowNeigh = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
		}
	}
	
	int noColRef = maxColRef-minColRef;
	int noRowRef = maxRowRef-minRowRef;
	
	int noColNeigh = maxColNeigh-minColNeigh;
	int noRowNeigh = maxRowNeigh-minRowNeigh;
	
	printf("matching img %d to %d\n", refIndex, neighIndex);
	float maxNoMatches = -1;
	
	//shift refIndex over neighIndex
	for(int i=-2/*-noColRef*/; i<=-1/*noColNeigh*/; i++)
	{
		for(int j=0/*-noRowRef*/; j<=0/*noRowNeigh*/; j++)
		{
			vector<KeypointMatch> tmpMatches;
			float noMatches = evaluateImgToImgGridAlignment(refIndex, neighIndex, templateId, refPlaneIndex, neighPlaneIndex, refGridIndex, neighGridIndex, i, j, tmpMatches);
			printf("%f matches for alignment (%d, %d)\n", noMatches, i, j);
			if(maxNoMatches == -1)
			{
				maxNoMatches = noMatches;
				col = i;
				row = j;
				bestMatches.clear();
				for(int m=0; m<tmpMatches.size(); m++)
				{
					bestMatches.push_back(tmpMatches[m]);
				}
			}
			else if(noMatches > maxNoMatches)
			{
				maxNoMatches = noMatches;
				col = i;
				row = j;
				bestMatches.clear();
				for(int m=0; m<tmpMatches.size(); m++)
				{
					bestMatches.push_back(tmpMatches[m]);
				}
			}
		}
	}
	
	//try without overlap between grids
	Matrix3f refM = rectifyingHomographies[refIndex][refPlaneIndex];
	Matrix3f refMInv = refM.getInverseMatrix();
	
	Matrix3f curM = rectifyingHomographies[neighIndex][neighPlaneIndex];
	Matrix3f curMInv = curM.getInverseMatrix();
	
	vector<KeypointMatch> tmpMatches = matches[refIndex][neighIndex - refIndex - 1];
	
	vector<Vec2f> aPts;
	vector<Vec2f> bPts;
	for(int m=0;m<tmpMatches.size(); m++)
	{
		Vec3f pos(keyInfo[refIndex][tmpMatches[m].m_idx1].x, keyInfo[refIndex][tmpMatches[m].m_idx1].y, 1.0);
		pos = refMInv*pos;
		pos[0] /= pos[2];
		pos[1] /= pos[2];
		aPts.push_back(shrink3To2(pos));
		
		Vec3f pos2(keyInfo[neighIndex][tmpMatches[m].m_idx2].x, keyInfo[neighIndex][tmpMatches[m].m_idx2].y, 1.0);
		pos2 = curMInv*pos2;
		pos2[0] /= pos2[2];
		pos2[1] /= pos2[2];
		bPts.push_back(shrink3To2(pos2));
		
	}
	vector<float> scale;
	vector<Vec2f> translation;
	vector<vector<int> > inliers;
	epGeo.estimateRectifiedTransformedRansacMatches(aPts.size(), aPts, bPts, 2048, 50, scale, translation, inliers);
	
	float threshold = 50.0;
	if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0] != 0.0 && 
	   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0]/2.0<threshold)
		threshold = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0]/2.0;
	if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1] != 0.0 && 
	   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1]/2.0<threshold)
		threshold = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1]/2.0;
	
	vector<vector<int> > tmpInliers;
	vector<int> newCol;
	vector<int> newRow;
	int bestCandidate;
	rf.translateTransformationToGridRANSAC(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex], 
										   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex], 
										   aPts, bPts, scale, translation, inliers, threshold, 
										   imageSizes[neighIndex][0], imageSizes[neighIndex][1], newCol, newRow, tmpInliers, bestCandidate);
	if(bestCandidate > -1)
	{
		printf("best:%d %d with %d\n", newCol[bestCandidate], newRow[bestCandidate], tmpInliers[bestCandidate].size());
	}
	//printf("bestCand:%d %d\n", newCol[bestCandidate], newRow[bestCandidate]);
	if(bestCandidate>-1 && (newCol[bestCandidate] < -noColRef || newCol[bestCandidate] > noColNeigh || newRow[bestCandidate] < -noRowRef || newRow[bestCandidate] > noRowNeigh))
	{
		if(tmpInliers[bestCandidate].size() > maxNoMatches)
		{
			bestMatches.clear();
			for(int m=0; m<tmpInliers[bestCandidate].size(); m++)
			{
				KeypointMatch kp = tmpMatches[tmpInliers[bestCandidate][m]];
				bestMatches.push_back(kp);
			}
			maxNoMatches = tmpInliers[bestCandidate].size();
			col = newCol[bestCandidate];
			row = newRow[bestCandidate];
		}
	}
	matchSize = maxNoMatches;
	printf("*********alignment from %d to %d is (%d, %d)\n", refIndex, neighIndex, col, row);
}

float BundlerManager::evaluateImgToImgGridAlignment(int refIndex, int neighIndex, int templateId, int refPlaneIndex, int neighPlaneIndex, int refGridIndex, int neighGridIndex, 
													int shiftCol, int shiftRow, vector<KeypointMatch> &tmpMatches)
{
	float totalMatches = 0;
	
    printf("a1\n");
	RepetitionFinder rf;
	EpipolarGeometry epGeo;
	
    printf("template:%d, refIndex:%d, refPlaneIndex:%d, refGridIndex:%d\n", templateId,refIndex,refPlaneIndex,refGridIndex);
    
	rf.fillMissingGridCellsTemporarily(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex], repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[0].size);
	int minCol, maxCol, minRow, maxRow;
	for(int c=0; c<repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells.size(); c++)
	{
		if(c==0)
		{
			minCol = maxCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			minRow = maxRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
		}
		else
		{
			if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber < minCol)
				minCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			else if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber > maxCol)
				maxCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			
			if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber < minRow)
				minRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
			else if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber > maxRow)
				maxRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
		}
	}
	
    printf("a2\n");
    printf("template:%d, neighIndex:%d, neighPlaneIndex:%d, neighGridIndex:%d\n", templateId,neighIndex,neighPlaneIndex,neighGridIndex);
	Matrix3f refM = rectifyingHomographies[refIndex][refPlaneIndex];
	Matrix3f refMInv = refM.getInverseMatrix();
	
	
	rf.fillMissingGridCellsTemporarily(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex], 
									   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[0].size);
	
    printf("a3\n");
    
	int minColCur, maxColCur, minRowCur, maxRowCur;
	for(int c=0; c<repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells.size(); c++)
	{
		if(c==0)
		{
			minColCur = maxColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			minRowCur = maxRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
		}
		else
		{
			if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber < minColCur)
				minColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			else if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber > maxColCur)
				maxColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			
			if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber < minRowCur)
				minRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
			else if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber > maxRowCur)
				maxRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
		}
	}
	
    printf("a4\n");
    
	Matrix3f curM = rectifyingHomographies[neighIndex][neighPlaneIndex];
	Matrix3f curMInv = curM.getInverseMatrix();
	
	minColCur = max(minCol+shiftCol, minColCur);
	maxColCur = min(maxCol+shiftCol, maxColCur);
	minRowCur = max(minRow+shiftRow, minRowCur);
	maxRowCur = min(maxRow+shiftRow, maxRowCur);
	
	Vec2f minCoord;
	Vec2f maxCoord;
	for(int c=0; c<repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells.size(); c++)
	{
		if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber >= minColCur && 
		   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber >= minRowCur
		   //repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber <= maxColCur && 
		   //repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber <= maxRowCur &&
		   //repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].temporary == false
		   )
		{
			Vec2f minC = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].center;
			minCoord = minC - Vec2f(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[0]/2.0, 
										repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[1]/2.0);
			maxCoord = minC + Vec2f(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[0]/2.0, 
											repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[1]/2.0);
			
		}
		if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber == maxColCur && 
		   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber == maxRowCur)
		{
			maxCoord = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].center;
			maxCoord = maxCoord + Vec2f(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[0]/2.0, 
										repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[1]/2.0);
		}
	}
	
    printf("a5\n");
    
	Vec2f minCoordRef;
	Vec2f maxCoordRef;
	for(int c=0; c<repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells.size(); c++)
	{
		if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber >= minColCur-shiftCol && 
		   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber >= minRowCur-shiftRow
		   //repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber <= maxColCur-shiftCol && 
		   //repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber <= maxRowCur-shiftRow &&
		   //repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].temporary == false
		   )
		{
			Vec2f minC = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].center;
			minCoordRef = minC - Vec2f(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[0]/2.0, 
											  repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[1]/2.0);
			maxCoordRef = minC + Vec2f(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[0]/2.0, 
											   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[1]/2.0);
		}
		if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber == maxColCur-shiftCol && 
		   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber == maxRowCur-shiftRow)
		{
			maxCoordRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].center;
			maxCoordRef = maxCoordRef + Vec2f(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[0]/2.0, 
											  repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[1]/2.0);
		}
	}
	
    printf("a6\n");
    
	//form new feature set
	vector<keypt_t> keyInfoCur;
	vector<keypt_t> keyInfoRef;
	vector<int> keyIdCur;
	vector<int> keyIdRef;
	vector<short> keysCur;
	vector<short> keysRef;
	
	//Img refImg, neighImg;
	//refImg.read(imageFileNames[refIndex]);
	//neighImg.read(imageFileNames[neighIndex]);
	
	//refImg.downScaleAndGaussianSmoothImage();
	//refImg.downScaleAndGaussianSmoothImage();
	//neighImg.downScaleAndGaussianSmoothImage();
	//neighImg.downScaleAndGaussianSmoothImage();
	
    printf("a7\n");
    
	for(int k=0; k<keyInfo[neighIndex].size(); k++)
	{
		Vec3f pos = Vec3f(keyInfo[neighIndex][k].x, keyInfo[neighIndex][k].y, 1.0);
		pos = curMInv*pos;
		pos[0] = pos[0]/pos[2];
		pos[1] = pos[1]/pos[2];
		bool found = false;
		if(pos[0]>minCoord[0] && pos[0] < maxCoord[0] && pos[1]>minCoord[1] && pos[1]<maxCoord[1])
		{
			found = true;
			continue;
		}
		
		//neighImg.setColor(keyInfo[neighIndex][k].x/4, keyInfo[neighIndex][k].y/4, Color(1.0, 0.0, 0.0));
		keyInfoCur.push_back(keyInfo[neighIndex][k]);
		keyIdCur.push_back(k);
		for(int d=0; d<128; d++)
		{
			keysCur.push_back(keys[neighIndex][k*128+d]);
		}
	}
    printf("a8\n");
    
	for(int k=0; k<keyInfo[refIndex].size(); k++)
	{
		Vec3f pos = Vec3f(keyInfo[refIndex][k].x, keyInfo[refIndex][k].y, 1.0);
		pos = refMInv*pos;
		pos[0] = pos[0]/pos[2];
		pos[1] = pos[1]/pos[2];
		bool found = false;
		if(pos[0]>minCoordRef[0] && pos[0] < maxCoordRef[0] && pos[1]>minCoordRef[1] && pos[1]<maxCoordRef[1])
		{
			found = true;
			continue;
		}
		
		//refImg.setColor(keyInfo[refIndex][k].x/4, keyInfo[refIndex][k].y/4, Color(1.0, 0.0, 0.0));
		keyInfoRef.push_back(keyInfo[refIndex][k]);
		keyIdRef.push_back(k);
		for(int d=0; d<128; d++)
		{
			keysRef.push_back(keys[refIndex][k*128+d]);
		}
	}
	
    printf("a9\n");
    
	//neighImg.write("/Users/ceylan/Desktop/neigh.png");
	//refImg.write("/Users/ceylan/Desktop/ref.png");
	
	vector<KeypointMatch> newMatches = kpMatcher->MatchKeys(keyInfoRef.size(), keysRef, keyInfoCur.size(), keysCur);
	
    printf("a10\n");
    
	vector<Vec2f> aPts;
	vector<Vec2f> bPts;
	for(int m=0;m<newMatches.size(); m++)
	{
		Vec3f pos(keyInfoRef[newMatches[m].m_idx1].x, keyInfoRef[newMatches[m].m_idx1].y, 1.0);
		pos = refMInv*pos;
		pos[0] /= pos[2];
		pos[1] /= pos[2];
		aPts.push_back(shrink3To2(pos));
		
		Vec3f pos2(keyInfoCur[newMatches[m].m_idx2].x, keyInfoCur[newMatches[m].m_idx2].y, 1.0);
		pos2 = curMInv*pos2;
		pos2[0] /= pos2[2];
		pos2[1] /= pos2[2];
		bPts.push_back(shrink3To2(pos2));
	}
	
	printf("a11\n");
    
	vector<float> scale;
	vector<Vec2f> translation;
	vector<vector<int> > inliers;
	epGeo.estimateRectifiedTransformedRansacMatches(aPts.size(), aPts, bPts, 2048, 50, scale, translation, inliers);
	
	float threshold = 50.0;
	if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0] != 0.0 && 
	   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0]/2.0<threshold)
		threshold = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0]/2.0;
	if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1] != 0.0 && 
	   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1]/2.0<threshold)
		threshold = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1]/2.0;
	
    printf("a12\n");
    
	vector<vector<int> > tmpInliers;
	vector<int> col;
	vector<int> row;
	int bestCandidate;
	rf.translateTransformationToGridRANSAC(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex], 
										   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex], 
										   aPts, bPts, scale, translation, inliers, threshold, 
										   imageSizes[neighIndex][0], imageSizes[neighIndex][1], col, row, tmpInliers, bestCandidate);
	
    printf("a13\n");
    
	printf("*******matching for %d and %d\n", shiftCol, shiftRow);
	for(int i=0; i<col.size(); i++)
	{
		printf("%d matches for %d, %d\n", tmpInliers[i].size(), col[i], row[i]);
		if(col[i] == shiftCol && row[i] == shiftRow)
		{
			for(int m=0; m<tmpInliers[i].size(); m++)
			{
				KeypointMatch kp;
				kp.m_idx1 = keyIdRef[newMatches[tmpInliers[i][m]].m_idx1];
				kp.m_idx2 = keyIdCur[newMatches[tmpInliers[i][m]].m_idx2];
				printf("%d to %d\n", kp.m_idx1, kp.m_idx2);
				tmpMatches.push_back(kp);
			}
		}
	}
	
	return tmpMatches.size();
}

float BundlerManager::findAlignmentByDiscardingGrids(int refIndex, int neighIndex, int templateId, int refPlaneIndex, int neighPlaneIndex, int refGridIndex, int neighGridIndex, 
													 int &shiftCol, int &shiftRow, vector<KeypointMatch> &tmpMatches)
{
	float totalMatches = 0;
	
	RepetitionFinder rf;
	EpipolarGeometry epGeo;
	
	rf.fillMissingGridCellsTemporarily(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex], repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[0].size);
	int minCol, maxCol, minRow, maxRow;
	for(int c=0; c<repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells.size(); c++)
	{
		if(c==0)
		{
			minCol = maxCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			minRow = maxRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
		}
		else
		{
			if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber < minCol)
				minCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			else if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber > maxCol)
				maxCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			
			if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber < minRow)
				minRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
			else if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber > maxRow)
				maxRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
		}
	}
	
	Matrix3f refM = rectifyingHomographies[refIndex][refPlaneIndex];
	Matrix3f refMInv = refM.getInverseMatrix();
	
	rf.fillMissingGridCellsTemporarily(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex], 
									   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[0].size);
	
	int minColCur, maxColCur, minRowCur, maxRowCur;
	for(int c=0; c<repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells.size(); c++)
	{
		if(c==0)
		{
			minColCur = maxColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			minRowCur = maxRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
		}
		else
		{
			if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber < minColCur)
				minColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			else if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber > maxColCur)
				maxColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			
			if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber < minRowCur)
				minRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
			else if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber > maxRowCur)
				maxRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
		}
	}
	
	Matrix3f curM = rectifyingHomographies[neighIndex][neighPlaneIndex];
	Matrix3f curMInv = curM.getInverseMatrix();
	
	Vec2f minCoord, maxCoord;
	for(int c=0; c<repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells.size(); c++)
	{
		if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber == minColCur && 
		   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber == minRowCur)
		{
			minCoord = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].center;
			minCoord = minCoord - Vec2f(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[0]/2.0, 
										repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[1]/2.0);
		}
		if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber == maxColCur && 
		   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber == maxRowCur)
		{
			maxCoord = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].center;
			maxCoord = maxCoord + Vec2f(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[0]/2.0, 
										repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[1]/2.0);
		}
	}
	
	Vec2f minCoordRef, maxCoordRef;
	for(int c=0; c<repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells.size(); c++)
	{
		if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber == minCol && 
		   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber == minRow)
		{
			minCoordRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].center;
			minCoordRef = minCoordRef - Vec2f(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[0]/2.0, 
											  repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[1]/2.0);
		}
		if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber == maxCol && 
		   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber == maxRow)
		{
			maxCoordRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].center;
			maxCoordRef = maxCoordRef + Vec2f(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[0]/2.0, 
											  repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[1]/2.0);
		}
	}
	
	//form new feature set
	vector<keypt_t> keyInfoCur;
	vector<keypt_t> keyInfoRef;
	vector<int> keyIdCur;
	vector<int> keyIdRef;
	vector<short> keysCur;
	vector<short> keysRef;
	
	for(int k=0; k<keyInfo[neighIndex].size(); k++)
	{
		Vec3f pos = Vec3f(keyInfo[neighIndex][k].x, keyInfo[neighIndex][k].y, 1.0);
		pos = curMInv*pos;
		pos[0] = pos[0]/pos[2];
		pos[1] = pos[1]/pos[2];
		if(pos[0]>minCoord[0] && pos[0] < maxCoord[0] && pos[1]>minCoord[1] && pos[1]<maxCoord[1])
		{
			continue;
		}
		
		keyInfoCur.push_back(keyInfo[neighIndex][k]);
		keyIdCur.push_back(k);
		for(int d=0; d<128; d++)
		{
			keysCur.push_back(keys[neighIndex][k*128+d]);
		}
	}
	for(int k=0; k<keyInfo[refIndex].size(); k++)
	{
		Vec3f pos = Vec3f(keyInfo[refIndex][k].x, keyInfo[refIndex][k].y, 1.0);
		pos = refMInv*pos;
		pos[0] = pos[0]/pos[2];
		pos[1] = pos[1]/pos[2];
		if(pos[0]>minCoordRef[0] && pos[0] < maxCoordRef[0] && pos[1]>minCoordRef[1] && pos[1]<maxCoordRef[1])
		{
			continue;
		}
		
		keyInfoRef.push_back(keyInfo[refIndex][k]);
		keyIdRef.push_back(k);
		for(int d=0; d<128; d++)
		{
			keysRef.push_back(keys[refIndex][k*128+d]);
		}
	}
	
	vector<KeypointMatch> newMatches = kpMatcher->MatchKeys(keyInfoRef.size(), keysRef, keyInfoCur.size(), keysCur);
	
	vector<Vec2f> aPts;
	vector<Vec2f> bPts;
	for(int m=0;m<newMatches.size(); m++)
	{
		Vec3f pos(keyInfoRef[newMatches[m].m_idx1].x, keyInfoRef[newMatches[m].m_idx1].y, 1.0);
		pos = refMInv*pos;
		pos[0] /= pos[2];
		pos[1] /= pos[2];
		aPts.push_back(shrink3To2(pos));
		
		Vec3f pos2(keyInfoCur[newMatches[m].m_idx2].x, keyInfoCur[newMatches[m].m_idx2].y, 1.0);
		pos2 = curMInv*pos2;
		pos2[0] /= pos2[2];
		pos2[1] /= pos2[2];
		bPts.push_back(shrink3To2(pos2));
		
	}
	vector<float> scale;
	vector<Vec2f> translation;
	vector<vector<int> > inliers;
	epGeo.estimateRectifiedTransformedRansacMatches(aPts.size(), aPts, bPts, 2048, 50, scale, translation, inliers);
	
	float threshold = 50.0;
	if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0] != 0.0 && 
	   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0]/2.0<threshold)
		threshold = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0]/2.0;
	if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1] != 0.0 && 
	   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1]/2.0<threshold)
		threshold = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1]/2.0;
	
	vector<vector<int> > tmpInliers;
	vector<int> col;
	vector<int> row;
	int bestCandidate;
	rf.translateTransformationToGridRANSAC(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex], 
										   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex], 
										   aPts, bPts, scale, translation, inliers, threshold, 
										   imageSizes[neighIndex][0], imageSizes[neighIndex][1], col, row, tmpInliers, bestCandidate);
	
	if(bestCandidate > -1)
	{
		shiftCol = col[bestCandidate];
		shiftRow = row[bestCandidate];
		for(int m=0; m<tmpInliers[bestCandidate].size(); m++)
		{
			KeypointMatch kp;
			kp.m_idx1 = keyIdRef[newMatches[tmpInliers[bestCandidate][m]].m_idx1];
			kp.m_idx2 = keyIdCur[newMatches[tmpInliers[bestCandidate][m]].m_idx2];
			tmpMatches.push_back(kp);
		}
		return tmpInliers[bestCandidate].size();
	}
	return 0;
}


float BundlerManager::findMatchesForAlignment(int refIndex, int neighIndex, int templateId, int refPlaneIndex, int neighPlaneIndex, int refGridIndex, int neighGridIndex, 
											  int shiftCol, int shiftRow, vector<KeypointMatch> &tmpMatches)
{
	float totalMatches = 0;
	
	RepetitionFinder rf;
	EpipolarGeometry epGeo;
	
	rf.fillMissingGridCellsTemporarily(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex], repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[0].size);
	int minCol, maxCol, minRow, maxRow;
	for(int c=0; c<repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells.size(); c++)
	{
		if(c==0)
		{
			minCol = maxCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			minRow = maxRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
		}
		else
		{
			if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber < minCol)
				minCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			else if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber > maxCol)
				maxCol = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber;
			
			if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber < minRow)
				minRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
			else if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber > maxRow)
				maxRow = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber;
		}
	}
	
	Matrix3f refM = rectifyingHomographies[refIndex][refPlaneIndex];
	Matrix3f refMInv = refM.getInverseMatrix();
	
	rf.fillMissingGridCellsTemporarily(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex], 
									   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[0].size);
	
	int minColCur, maxColCur, minRowCur, maxRowCur;
	for(int c=0; c<repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells.size(); c++)
	{
		if(c==0)
		{
			minColCur = maxColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			minRowCur = maxRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
		}
		else
		{
			if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber < minColCur)
				minColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			else if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber > maxColCur)
				maxColCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber;
			
			if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber < minRowCur)
				minRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
			else if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber > maxRowCur)
				maxRowCur = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber;
		}
	}
	
	Matrix3f curM = rectifyingHomographies[neighIndex][neighPlaneIndex];
	Matrix3f curMInv = curM.getInverseMatrix();
	
	Vec2f minCoord, maxCoord;
	for(int c=0; c<repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells.size(); c++)
	{
		if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber == minColCur && 
		   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber == minRowCur)
		{
			minCoord = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].center;
			minCoord = minCoord - Vec2f(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[0]/2.0, 
										repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[1]/2.0);
		}
		if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].colNumber == maxColCur && 
		   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].rowNumber == maxRowCur)
		{
			maxCoord = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].center;
			maxCoord = maxCoord + Vec2f(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[0]/2.0, 
										repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].gridCells[c].size[1]/2.0);
		}
	}
	
	Vec2f minCoordRef, maxCoordRef;
	for(int c=0; c<repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells.size(); c++)
	{
		if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber == minCol && 
		   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber == minRow)
		{
			minCoordRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].center;
			minCoordRef = minCoordRef - Vec2f(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[0]/2.0, 
											  repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[1]/2.0);
		}
		if(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].colNumber == maxCol && 
		   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].rowNumber == maxRow)
		{
			maxCoordRef = repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].center;
			maxCoordRef = maxCoordRef + Vec2f(repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[0]/2.0, 
											  repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex].gridCells[c].size[1]/2.0);
		}
	}
	
	//form new feature set
	vector<keypt_t> keyInfoCur;
	vector<keypt_t> keyInfoRef;
	vector<int> keyIdCur;
	vector<int> keyIdRef;
	vector<short> keysCur;
	vector<short> keysRef;
	
	for(int k=0; k<keyInfo[neighIndex].size(); k++)
	{
		Vec3f pos = Vec3f(keyInfo[neighIndex][k].x, keyInfo[neighIndex][k].y, 1.0);
		pos = curMInv*pos;
		pos[0] = pos[0]/pos[2];
		pos[1] = pos[1]/pos[2];
		keyInfoCur.push_back(keyInfo[neighIndex][k]);
		keyIdCur.push_back(k);
		for(int d=0; d<128; d++)
		{
			keysCur.push_back(keys[neighIndex][k*128+d]);
		}
	}
	for(int k=0; k<keyInfo[refIndex].size(); k++)
	{
		Vec3f pos = Vec3f(keyInfo[refIndex][k].x, keyInfo[refIndex][k].y, 1.0);
		pos = refMInv*pos;
		pos[0] = pos[0]/pos[2];
		pos[1] = pos[1]/pos[2];
		
		keyInfoRef.push_back(keyInfo[refIndex][k]);
		keyIdRef.push_back(k);
		for(int d=0; d<128; d++)
		{
			keysRef.push_back(keys[refIndex][k*128+d]);
		}
	}
	
	vector<KeypointMatch> newMatches = kpMatcher->MatchKeys(keyInfoRef.size(), keysRef, keyInfoCur.size(), keysCur);
	
	vector<Vec2f> aPts;
	vector<Vec2f> bPts;
	for(int m=0;m<newMatches.size(); m++)
	{
		Vec3f pos(keyInfoRef[newMatches[m].m_idx1].x, keyInfoRef[newMatches[m].m_idx1].y, 1.0);
		pos = refMInv*pos;
		pos[0] /= pos[2];
		pos[1] /= pos[2];
		aPts.push_back(shrink3To2(pos));
		
		Vec3f pos2(keyInfoCur[newMatches[m].m_idx2].x, keyInfoCur[newMatches[m].m_idx2].y, 1.0);
		pos2 = curMInv*pos2;
		pos2[0] /= pos2[2];
		pos2[1] /= pos2[2];
		bPts.push_back(shrink3To2(pos2));
		
	}
	vector<float> scale;
	vector<Vec2f> translation;
	vector<vector<int> > inliers;
	epGeo.estimateRectifiedTransformedRansacMatches(aPts.size(), aPts, bPts, 2048, 50, scale, translation, inliers);
	
	float threshold = 50.0;
	if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0] != 0.0 && 
	   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0]/2.0<threshold)
		threshold = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].xTrans[0]/2.0;
	if(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1] != 0.0 && 
	   repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1]/2.0<threshold)
		threshold = repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex].yTrans[1]/2.0;
	
	vector<vector<int> > tmpInliers;
	vector<int> col;
	vector<int> row;
	int bestCandidate;
	rf.translateTransformationToGridRANSAC(repeatingGrids[templateId][neighIndex][neighPlaneIndex][neighGridIndex], 
										   repeatingGrids[templateId][refIndex][refPlaneIndex][refGridIndex], 
										   aPts, bPts, scale, translation, inliers, threshold, 
										   imageSizes[neighIndex][0], imageSizes[neighIndex][1], col, row, tmpInliers, bestCandidate);
	
	for(int i=0; i<col.size(); i++)
	{
		if(col[i] == shiftCol && row[i] == shiftRow)
		{
			for(int m=0; m<tmpInliers[i].size(); m++)
			{
				KeypointMatch kp;
				kp.m_idx1 = keyIdRef[newMatches[tmpInliers[i][m]].m_idx1];
				kp.m_idx2 = keyIdCur[newMatches[tmpInliers[i][m]].m_idx2];
				tmpMatches.push_back(kp);
			}
			return tmpInliers[bestCandidate].size();
		}
	}
	return 0;
}

void BundlerManager::computeTracks()
{
	tracks.clear();
	if(kpMatcher == NULL)
		kpMatcher = new KeypointMatcher();
	kpMatcher->computeTracks(keyInfo, matches, tracks);	
}

int BundlerManager::getNoOfTrack()
{
	return tracks.size();
}

int BundlerManager::getNoOfCorrespondencesInTrack(int trackIndex)
{
	if(trackIndex >= tracks.size())
		return 0;
	else
		tracks[trackIndex].correspondences.size();
}

Correspondence BundlerManager::getCorrespondenceInTrack(int trackIndex, int imageIndex)
{
	return tracks[trackIndex].correspondences[imageIndex];
		
}

void BundlerManager::formMatchesFromTracks()
{
	BundlerOutputParser *bundleParser = new BundlerOutputParser();
	
	string fname = "/Users/ceylan/MVSReconstruction/data/MVSDataSets/UCLPubBuilding/fullSet/FixedFocal/bundle.out";
	
	bundleParser->readBundlerOutputFile(fname.c_str(), 2160, 2880);
	
	matches.clear();
	matches.resize(imageFileNames.size());
	fundMatrices.clear();
	fundMatrices.resize(imageFileNames.size());
	
	for(int i=0; i<imageFileNames.size(); i++)
	{
		matches[i].resize(imageFileNames.size()-i-1);
		fundMatrices[i].resize(imageFileNames.size()-i-1);
	}
	
	int totalNoTracks = bundleParser->getNoOfTrack();
	for(int i=0; i<totalNoTracks; i++)
	{
		int visibleImageNo = bundleParser->getNoOfCorrespondencesInTrack(i);
		
		vector<Correspondence> c;
		for(int j=0; j<visibleImageNo; j++)
		{
			c.push_back(bundleParser->getCorrespondenceInTrack(i, j));
		}
		
		for(int j=0; j<c.size(); j++)
		{
			int imageIndex1 = c[j].imageIndex;
			int keyIndex1 = c[j].keyIndex;
			
			for(int k=j+1; k<c.size(); k++)
			{
				int imageIndex2 = c[k].imageIndex;
				int keyIndex2 = c[k].keyIndex;
				
				if(imageIndex1 < imageIndex2)
				{
					KeypointMatch kp;
					kp.m_idx1 = keyIndex1;
					kp.m_idx2 = keyIndex2;
					matches[imageIndex1][imageIndex2-imageIndex1-1].push_back(kp);
				}
				else
				{
					KeypointMatch kp;
					kp.m_idx1 = keyIndex2;
					kp.m_idx2 = keyIndex1;
					matches[imageIndex2][imageIndex1-imageIndex2-1].push_back(kp);
				}
			}
		}
	}
	
	EpipolarGeometry eg;
	string fundFile = mainDirectory + "/fundMatrix.txt";
	ifstream fin(fundFile.c_str(), ios::in);
	if(fin)
	{
		int noMatrices;
		fin >> noMatrices;
		
		for(int i=0; i<imageFileNames.size(); i++)
		{
			for(int j=i+1; j<imageFileNames.size(); j++)
			{
				Matrix3f F;
				F.setIdentity();
				
				fin >> F[0][0] >> F[1][0] >> F[2][0];
				fin >> F[0][1] >> F[1][1] >> F[2][1];
				fin >> F[0][2] >> F[1][2] >> F[2][2];
				
				fundMatrices[i][j-i-1] = F;
			}
		}
	}
	else
	{
		ofstream fout(fundFile.c_str(), ios::out);
		fout << (imageFileNames.size()*(imageFileNames.size()-1)/2) << endl;
		
		for(int i=0; i<imageFileNames.size(); i++)
		{
			for(int j=i+1; j<imageFileNames.size(); j++)
			{
				printf("matching image %d to image %d\n", i, j);
				Matrix3f F;
				F.setIdentity();
				if(matches[i][j-i-1].size() >= 8)
					eg.estimateFMatrix(keyInfo[i], keyInfo[j], matches[i][j-i-1], 1024, 20.0, F);
				fout << F[0][0] << " " << F[1][0] << " " << F[2][0] << endl;
				fout << F[0][1] << " " << F[1][1] << " " << F[2][1] << endl;
				fout << F[0][2] << " " << F[1][2] << " " << F[2][2] << endl;
				fundMatrices[i][j-i-1] = F;
			}
		}
		fout.close();
	}
	
}

void BundlerManager::runBundleAdjustment()
{
	if(kpMatcher == NULL)
		kpMatcher = new KeypointMatcher();
	
	tracks.clear();
	//kpMatcher->computeTracks(keyInfo, matches, tracks);
	
	if(geoManager == NULL)
		geoManager = new GeometryManager();
	
	bool removeBadMatches = true;
	
	for(int i=0; i<transformations.size(); i++)
		transformations[i].clear();
	transformations.clear();
	
	int numImages = imageFileNames.size();
	
	for(int i=0; i<numImages; i++)
	{
		transformations.push_back(std::vector<Transformation>(numImages));
	}
	
	//geoManager->computeEpipolarGeometry(transformations, removeBadMatches, keyInfo, matches);
	
	//return;
	
	string listFile = mainDirectory + "/list.txt";
	FILE *fl = fopen(listFile.c_str(), "w");
	for(int i=0; i<imageFileNames.size(); i++)
	{
		fprintf(fl, "%s\n", imageFileNames[i].c_str());
	}
	fclose(fl);
	
	string filename = mainDirectory + "/options.txt";
	printf("file:%s\n", filename.c_str());
	
	FILE *f = fopen(filename.c_str(), "w");
	
	string keyDirectory = mainDirectory + "/matches";
	fprintf(f, "%s%s\n",  "--key_dir ", keyDirectory.c_str());
	
	string matchFile = mainDirectory + "/matches_init.txt";
	fprintf(f, "%s%s\n",  "--match_table ", matchFile.c_str());
	
	string bundleDir = mainDirectory + "/bundle";
	fprintf(f, "%s\n", "--output bundle.out");
	fprintf(f, "%s%s\n", "--output_dir ", bundleDir.c_str());
	bundlerOutputPath = bundleDir + "/bundle.out";
	
	fprintf(f, "%s\n", "--fixed_focal_length");
	//fprintf(f, "%s\n", "--variable_focal_length");
	//fprintf(f, "%s\n", "--use_focal_estimate");
	//fprintf(f, "%s\n", "--trust_focal_estimate");
	//fprintf(f, "%s\n", "--constrain_focal");
	//fprintf(f, "%s\n", "--constrain_focal_weight 0.1");
	//fprintf(f, "%s\n", "--initial_focal_length 1493.0");
	//fprintf(f, "%s\n", "--estimate_distortion");
	//fprintf(f, "%s %s\n", "--intrinsics", std::string(mainDirectory + "/intrinsics.txt").c_str());
	fprintf(f, "%s\n", "--run_bundle");
	
	fclose(f);
	
	int argc = 3;
	char **argv = new char*[3];
	argv[0] = (char*) listFile.c_str();
	argv[1] = "--options_file";
	argv[2] = (char*)filename.c_str();
	
	BundlerApp *bundler_app = new BundlerApp();
    bundler_app->argc = argc;
    bundler_app->argv = argv;
	
	//check if there's a repetition file
	string gridFile = mainDirectory + "/matches/repetitions.txt";
	ifstream fin(gridFile.c_str());
	if(fin)
	{
		int noImg, noGrids;
		int *noRow, *noCol;
		int found;
		vector<vector<vector<int> > > gridMask;
		fin >> noImg >> noGrids;
		noRow = new int[noGrids];
		noCol = new int[noGrids];
		gridMask.resize(noGrids);
		for(int g=0; g<noGrids; g++)
		{
			fin >> noRow[g] >> noCol[g];
			gridMask[g].resize(noImg);
		
			for(int i=0; i<noImg; i++)
			{
				for(int r=0; r<noRow[g]; r++)
				{
					for(int c=0; c<noCol[g]; c++)
					{
						fin >> found;
						gridMask[g][i].push_back(found);
					}
				}
			}
		}
		
		fin.close();
		
		vector<Vec3f> refPoint, horTrans, verTrans;
		
		bundler_app->OnInitWithGrid(noGrids, noRow, noCol, gridMask);
		
		string gridFile = mainDirectory + "/matches/grid.txt";
		ofstream gridFin(gridFile.c_str(), ios::out);
		gridFin << noGrids << endl;
		float x,y,z;
		for(int g=0; g<noGrids; g++)
		{
			bundler_app->getGridRefPoint(g, x, y, z);
			refPoint.push_back(Vec3f(x,y,z));
			gridFin << x << " " << y << " " << z << " ";
			bundler_app->getGridHorTrans(g, x, y, z);
			horTrans.push_back(Vec3f(x,y,z));
			gridFin << x << " " << y << " " << z << " ";
			bundler_app->getGridVerTrans(g, x, y, z);
			verTrans.push_back(Vec3f(x,y,z));
			gridFin << x << " " << y << " " << z << endl;
		}
		gridFin.close();
	}
	else
	{
		bundler_app->OnInit();
		printf("ERROR: Repetition mask file not found.\n");
		return ;
	}
	
	points3D.clear();
	
	camMatrices.clear();
	camMatrices.resize(numImages);
	
	for(int i=0; i<pointProjections.size(); i++)
		pointProjections[i].clear();
	for(int i=0; i<pointProjectionColors.size(); i++)
		pointProjectionColors[i].clear();
	
	pointProjections.clear();
	pointProjections.resize(numImages);
	pointProjectionColors.clear();
	pointProjectionColors.resize(numImages);
	
	int no3DPoints = bundler_app->m_point_data.size();
	
	vector<vector<int> > corrFeatures;
	
	for(int i=0; i<no3DPoints; i++)
	{
		points3D.push_back(Vec3f((bundler_app->m_point_data)[i].m_pos[0],(bundler_app->m_point_data)[i].m_pos[1],-(bundler_app->m_point_data)[i].m_pos[2]));
		
		printf("******point %d*******\n", i);
		vector<int> visFeatures;
		visFeatures.resize(numImages, -1);
		
		for(int v=0; v<(bundler_app->m_point_data)[i].m_views.size(); v++)
		{
			int keyNo = (bundler_app->m_point_data)[i].m_views[v].second;
			int imgNo = (bundler_app->m_point_data)[i].m_views[v].first;
			visFeatures[imgNo] = keyNo;
			printf("visible in img %d as feature %d\n", imgNo, keyNo);
		}
		
		corrFeatures.push_back(visFeatures);
		
	}
	
	for(int i=0; i<numImages; i++)
	{
		double focal = bundler_app->m_image_data[i].m_camera.m_focal;
		double width = bundler_app->m_image_data[i].GetWidth();
		double height = bundler_app->m_image_data[i].GetHeight();
		double *R = bundler_app->m_image_data[i].m_camera.m_R;
		double *t = bundler_app->m_image_data[i].m_camera.m_t;
		
		Matrix4f K;
		K.setIdentity();
		K[0][0] = -focal; K[2][0] = 0.5*width-0.5;
		K[1][1] = focal; K[2][1] = 0.5*height-0.5;
		
		Matrix4f P;
		P.setIdentity();
		P[0][0] = R[0]; P[1][0] = R[1]; P[2][0] = -R[2]; P[3][0] = t[0];
		P[0][1] = R[3]; P[1][1] = R[4]; P[2][1] = -R[5]; P[3][1] = t[1];
		P[0][2] = -R[6]; P[1][2] = -R[7]; P[2][2] = R[8]; P[3][2] = -t[2];
		
		P = K*P;
		P = P*-1.0;
		camMatrices[i] = P;
		
		//printf("focal:%f\n", focal);
		//printf("R:%f %f %f %f %f %f %f %f %f\n", R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
		//printf("t:%f %f %f\n", t[0], t[1], t[2]);
		//printf("PMat:\n");
		//printf("%f %f %f %f\n", P[0][0], P[1][0], P[2][0], P[3][0]);
		//printf("%f %f %f %f\n", P[0][1], P[1][1], P[2][1], P[3][1]);
		//printf("%f %f %f %f\n", P[0][2], P[1][2], P[2][2], P[3][2]);
		
		for(int j=0; j<points3D.size(); j++)
		{
			//printf("point:%f %f %f\n", points3D[j][0], points3D[j][1], points3D[j][2]);
			Vec3f pos = shrink4To3(P*expand3To4(points3D[j]));
			Vec3uc c;
			//if(corrFeatures[j][i] != -1.0)
			{
				printf("imageNo:%d, featureNo: %d\n", i, corrFeatures[j][i]);
				
				//Vec2f featurePos = Vec2f(keyInfo[i][corrFeatures[j][i]].x, keyInfo[i][corrFeatures[j][i]].y);
				//float dist = (Vec2f(pos[0]/pos[2], pos[1]/pos[2])-featurePos).length();
				//Vec3uc c = MathUtils::generateColorFromValue(dist, 0, 50.0);
				Vec3uc c(0,0,255);
				//printf(" dist:%f\n", dist);
				
				if(pos[2] != 0.0)
				{
					pos[0] /= pos[2]; pos[1] /= pos[2];
				
				
					if(pos[0]>=0.0 && pos[0]<width && pos[1]>=0.0 && pos[1]<height)
					{
						pointProjections[i].push_back(Vec2f(pos[0], pos[1]));
						pointProjectionColors[i].push_back(Vec3f(((float)(c[0]))/255.0, ((float)(c[1]))/255.0, ((float)(c[2]))/255.0));
						//printf("proj:%f %f\n", pos[0], pos[1]);
					}
				}
			}
		}
	}
}

void BundlerManager::writeCalibrationResults()
{
	if(boParser == NULL)
	{
		boParser = new BundlerOutputParser();
	}
	//int width = imageSizes[0][0];
	//int height = imageSizes[0][1];
	//boParser->readBundlerOutputFile(bundlerOutputPath.c_str(), width, height);
	int width = 1000;
    int height = 1000;
    boParser->readBundlerOutputFile("/Users/ceylan/MVSReconstruction/data/MVSDataSets/Roman_Forum/output/gt_bundle.out", width, height);
	string calibResultPath = mainDirectory + "/bundler";
	cameras.clear();
	boParser->writeCamAndImage(calibResultPath.c_str(), imageFileNames, cameras);
	
	generateCameraFrustums();
}

void BundlerManager::readCalibrationResults()
{
	if(boParser == NULL)
	{
		boParser = new BundlerOutputParser();
	}
	
	string calibResultPath = mainDirectory + "/pmvs";
	cameras.clear();
	boParser->readCalibration(calibResultPath.c_str(), imageFileNames.size(), cameras);
	
	generateCameraFrustums();
}


void BundlerManager::generateCameraFrustums()
{
	Mesh3D camMesh;
	
	//int camOrder[] = {17, 20, 18, 19, 21, 22, 23, 9, 8, 13, 6, 5, 11, 4, 7, 12, 3, 2, 1, 0, 10, 25, 14, 24, 15, 16};
	
	for(int i=0; i<cameras.size(); i++)
	{
		int index = i;
		
		//int index = camOrder[i];
		
		Vec3f center = shrink4To3(cameras[index].getCenter());
        
        printf("center:%f %f %f\n", center[0], center[1], center[2]);
        
		Vec3f ll = cameras[index].unprojectPixelwithDepth(Vec3f(0.0, imageSizes[index][1], 1), 0);
		Vec3f lr = cameras[index].unprojectPixelwithDepth(Vec3f(imageSizes[index][0], imageSizes[index][1], 1), 0);
		Vec3f ul = cameras[index].unprojectPixelwithDepth(Vec3f(0.0, 0.0, 1), 0);
		Vec3f ur = cameras[index].unprojectPixelwithDepth(Vec3f(imageSizes[index][0], 0.0, 1), 0);

		Vec3f dir;
		dir = ll - center;
		dir.normalize();
		ll = center + dir*0.4;
		
		dir = lr - center;
		dir.normalize();
		lr = center + dir*0.4;
		
		dir = ul - center;
		dir.normalize();
		ul = center + dir*0.4;
		
		dir = ur - center;
		dir.normalize();
		ur = center + dir*0.4;
		
		camMesh.addPoint(center);
		camMesh.addPoint(ll);
		camMesh.addPoint(lr);
		camMesh.addPoint(ul);
		camMesh.addPoint(ur);
		
		camMesh.addFace(i*5, i*5+3, i*5+1);
		camMesh.addFace(i*5, i*5+1, i*5+2);
		camMesh.addFace(i*5, i*5+2, i*5+4);
		camMesh.addFace(i*5, i*5+4, i*5+3);
		
		camMesh.addFace(i*5+1, i*5+3, i*5+4);
		camMesh.addFace(i*5+4, i*5+2, i*5+1);
		
		stringstream ss;
		ss << mainDirectory << "/pmvs/cameras" << i << ".obj";
		
		string filename;
		ss >> filename;
		
		camMesh.writeMesh(filename.c_str());
	}
	
}

Matrix3f getRotationMatrix(float angle, int axis)
{
	Matrix3f rot;
	
	angle = angle * M_PI / 180.0;
	
	switch(axis)
	{
		case 1:
			//x axis
			rot[0][0] = 1;
			rot[1][0] = 0;
			rot[2][0] = 0;
			rot[0][1] = 0; 
			rot[1][1] = cos(angle);
			rot[2][1] = sin(angle); 
			rot[0][2] = 0;
			rot[1][2] = -sin(angle);
			rot[2][2] = cos(angle);
			break;
			
		case 2:
			//y axis
			rot[0][0] = cos(angle);
			rot[1][0] = 0;
			rot[2][0] = -sin(angle);
			rot[0][1] = 0; 
			rot[1][1] = 1;
			rot[2][1] = 0; 
			rot[0][2] = sin(angle);
			rot[1][2] = 0;
			rot[2][2] = cos(angle);
			break;
			
		case 3:
			//z axis
			rot[0][0] = cos(angle);
			rot[1][0] = sin(angle);
			rot[2][0] = 0;
			rot[0][1] = -sin(angle); 
			rot[1][1] = cos(angle);
			rot[2][1] = 0; 
			rot[0][2] = 0;
			rot[1][2] = 0;
			rot[2][2] = 1;
			break;
			
		default:
			rot.setIdentity();
	}
	return rot;
}

void BundlerManager::bringToSameCoordinateFrame(vector<Vec2f> imageCoordinates1, vector<Vec2f> imageCoordinates2, int cam1, int cam2)
{
	EpipolarGeometry eg;
	
	/*Matrix3f K;
	K.setIdentity();
	K[0][0] = 2209.091;
	K[1][1] = 2209.091;
	K[2][0] = 1080;
	K[2][1] = 810;
	
	Matrix3f Kinv = K.getInverseMatrix();
	
	Matrix3f rotX180;
	rotX180[0][0] = 1;
	rotX180[1][0] = 0;
	rotX180[2][0] = 0;
	rotX180[0][1] = 0; 
	rotX180[1][1] = cos(M_PI);
	rotX180[2][1] = sin(M_PI); 
	rotX180[0][2] = 0;
	rotX180[1][2] = -sin(M_PI);
	rotX180[2][2] = cos(M_PI);
	
	float apertureX = 1.417 * 25.4;
	float apertureY = 0.945 * 25.4;

	for(int i=0; i<cameras.size(); i++)
	{
		Matrix4f myProj = expand3To4(Kinv) * cameras[i].getProjectionMatrix(0);
		Matrix3f rot = shrink4To3(myProj);
		Vec3f t = shrink4To3(myProj[3]);
		
		Matrix3f rotCheck = K * rot;
		Vec3f tCheck = K*t;
		
		printf("check:\n");
		printf("%f %f %f %f\n", rotCheck[0][0], rotCheck[1][0], rotCheck[2][0], tCheck[0]);
		printf("%f %f %f %f\n", rotCheck[0][1], rotCheck[1][1], rotCheck[2][1], tCheck[1]);
		printf("%f %f %f %f\n", rotCheck[0][2], rotCheck[1][2], rotCheck[2][2], tCheck[2]);
		
		Matrix3f mayaRot = rotX180 * rot;
		Vec3f mayaCenter = shrink4To3(cameras[i].getCenter());
		
		float f = 2209.091 * apertureY / 1620.0;
		
		float angleY1 = -asin(mayaRot[0][2]);
		float angleY2 = M_PI - angleY1;
		float angleX1 = atan2(mayaRot[1][2]/cos(angleY1), mayaRot[2][2]/cos(angleY1));
		float angleX2 = atan2(mayaRot[1][2]/cos(angleY2), mayaRot[2][2]/cos(angleY2));
		float angleZ1 = atan2(mayaRot[0][1]/cos(angleY1), mayaRot[0][0]/cos(angleY1));
		float angleZ2 = atan2(mayaRot[0][1]/cos(angleY2), mayaRot[0][0]/cos(angleY2));
		
		mayaRot = getRotationMatrix(3, angleZ1) * getRotationMatrix(2, angleY1) * getRotationMatrix(1, angleX1);
		
		angleX1 = angleX1 * 180.0 / M_PI;
		angleY1 = angleY1 * 180.0 / M_PI;
		angleZ1 = angleZ1 * 180.0 / M_PI;
		
		angleX2 = angleX2 * 180.0 / M_PI;
		angleY2 = angleY2 * 180.0 / M_PI;
		angleZ2 = angleZ2 * 180.0 / M_PI;
		
		printf("angleX:%f, angleY:%f, angleZ:%f\n", angleX1, angleY1, angleZ1);
		printf("angleX:%f, angleY:%f, angleZ:%f\n", angleX2, angleY2, angleZ2);
		printf("cam:%f %f %f\n", mayaCenter[0], mayaCenter[1], mayaCenter[2]);
		printf("focal:%f\n", f);
	}
	
	return ;*/
	
    vector<double> focalLengths;
    focalLengths.push_back(1.3128579625e+003);
    focalLengths.push_back(1.3158687288e+003);
    focalLengths.push_back(1.3273590562e+003);
    focalLengths.push_back(1.3235880180e+003);
    focalLengths.push_back(1.3219821470e+003);
    focalLengths.push_back(1.3124006149e+003);
    focalLengths.push_back(1.3002940992e+003);
    focalLengths.push_back(1.3000039682e+003);
    focalLengths.push_back(1.2879774932e+003);
    focalLengths.push_back(1.3067531289e+003);
    focalLengths.push_back(1.2930380430e+003);
    focalLengths.push_back(1.2905685892e+003);
    focalLengths.push_back(1.2936679070e+003);
    focalLengths.push_back(1.3033424923e+003);
    focalLengths.push_back(1.2905611885e+003);
    focalLengths.push_back(1.2951620419e+003);
    focalLengths.push_back(1.2894999531e+003);
    focalLengths.push_back(1.2840406069e+003);
    focalLengths.push_back(1.2760146531e+003);
    focalLengths.push_back(1.2950846260e+003);
    focalLengths.push_back(1.2920206510e+003);
    focalLengths.push_back(1.2832281717e+003);
    focalLengths.push_back(1.2956278166e+003);
    focalLengths.push_back(1.2769583137e+003);
    
    /*focalLengths.push_back(1.1240600506e+003);
    focalLengths.push_back(1.1201433505e+003);
    focalLengths.push_back(1.1276591295e+003);
    focalLengths.push_back(1.1194707183e+003);
    focalLengths.push_back(1.1165789221e+003);
    focalLengths.push_back(1.1187540691e+003);
    focalLengths.push_back(1.1170930777e+003);
    focalLengths.push_back(1.1138550662e+003);
    focalLengths.push_back(1.1081082617e+003);
    focalLengths.push_back(1.1114630122e+003);
    focalLengths.push_back(1.1083711686e+003);
    focalLengths.push_back(1.1105822063e+003);
    focalLengths.push_back(1.1094631081e+003);
    focalLengths.push_back(1.1059183918e+003);
    focalLengths.push_back(1.1070295356e+003);
    focalLengths.push_back(1.1230537894e+003);
    focalLengths.push_back(1.1247731578e+003);
    focalLengths.push_back(1.1372072130e+003);
    focalLengths.push_back(1.1373520304e+003);
    focalLengths.push_back(1.1457767858e+003);
    focalLengths.push_back(1.1490128562e+003);
    focalLengths.push_back(1.1327160735e+003);
    focalLengths.push_back(1.0970939228e+003);
    focalLengths.push_back(1.0907023635e+003);
    focalLengths.push_back(1.0829035528e+003);
    focalLengths.push_back(1.0881500438e+003);
    focalLengths.push_back(1.0721658243e+003);
    */
    
    vector<Matrix3f> KMats;
    vector<Matrix3f> KInvMats;
    for(int i=0; i<focalLengths.size(); i++)
    {
        Matrix3f K;
        K.setIdentity();
        
        if(focalLengths[i] != 0.0)
        {
            K[0][0] = focalLengths[i];
            K[1][1] = focalLengths[i];
            K[2][0] = 1067/2.0;
            K[2][1] = 800/2.0;
        }
        
        Matrix3f Kinv = K.getInverseMatrix();
        
        KMats.push_back(K);
        KInvMats.push_back(Kinv);
    }
    
	//current camera coordinate system
	//line1
	Vec3f s1cur = eg.computeClosest3DPoint(imageCoordinates1[0], imageCoordinates2[0], cameras[cam1], cameras[cam2], 0);
	Vec3f e1cur = eg.computeClosest3DPoint(imageCoordinates1[1], imageCoordinates2[1], cameras[cam1], cameras[cam2], 0);
		
	//line2
	Vec3f s2cur = eg.computeClosest3DPoint(imageCoordinates1[2], imageCoordinates2[2], cameras[cam1], cameras[cam2], 0);
	Vec3f e2cur = eg.computeClosest3DPoint(imageCoordinates1[3], imageCoordinates2[3], cameras[cam1], cameras[cam2], 0);
	
	vector<Vec3f> refPoints;
	refPoints.push_back(s1cur);
	refPoints.push_back(e1cur);
	refPoints.push_back(s2cur);
	refPoints.push_back(e2cur);
	
	string camDirectory = mainDirectory + "/comparisonJiang2012";
	if(boParser == NULL)
	{
		boParser = new BundlerOutputParser();
	}
	vector<PerspectiveCamera> newCameras;
	boParser->readCalibration(camDirectory.c_str(), imageFileNames.size(), newCameras);
	
	
	//line1
	Vec3f s1 = eg.computeClosest3DPoint(imageCoordinates1[0]/2.7, imageCoordinates2[0]/2.7, newCameras[cam1], newCameras[cam2], 0);
	Vec3f e1 = eg.computeClosest3DPoint(imageCoordinates1[1]/2.7, imageCoordinates2[1]/2.7, newCameras[cam1], newCameras[cam2], 0);
	
	//line2
	Vec3f s2 = eg.computeClosest3DPoint(imageCoordinates1[2]/2.7, imageCoordinates2[2]/2.7, newCameras[cam1], newCameras[cam2], 0);
	Vec3f e2 = eg.computeClosest3DPoint(imageCoordinates1[3]/2.7, imageCoordinates2[3]/2.7, newCameras[cam1], newCameras[cam2], 0);
		
	
	float scale = (s1cur-e1cur).length() / (s1 - e1).length();
		
	printf("scale:%f\n", scale);
		
	PointCloud *pc = new PointCloud("/Users/ceylan/MVSReconstruction/data/MVSDataSets/ucl5/comparisonJiang2012/dense.ply");
	pc->scaleData(Vec3f(scale, scale, scale));
	
	pc->writeMesh("/Users/ceylan/Desktop/deneme.ply");
	
	//apply scaling
	Vec3f newCenter;
	Matrix4f scaleMat;
		
	newCenter = shrink4To3(newCameras[cam1].getCenter()) * scale;
	newCenter = -shrink4To3(newCameras[cam1].getProjectionMatrix(0)) * newCenter;
	scaleMat = newCameras[cam1].getProjectionMatrix(0);
	scaleMat[3][0] = newCenter[0]; scaleMat[3][1] = newCenter[1]; scaleMat[3][2] = newCenter[2];
	newCameras[cam1].setProjectionMatrix(scaleMat);
		
	newCenter = shrink4To3(newCameras[cam2].getCenter()) * scale;
	newCenter = -shrink4To3(newCameras[cam2].getProjectionMatrix(0)) * newCenter;
	scaleMat = newCameras[cam2].getProjectionMatrix(0);
	scaleMat[3][0] = newCenter[0]; scaleMat[3][1] = newCenter[1]; scaleMat[3][2] = newCenter[2];
	newCameras[cam2].setProjectionMatrix(scaleMat);
		
	//line1
	s1 = eg.computeClosest3DPoint(imageCoordinates1[0]/2.7, imageCoordinates2[0]/2.7, newCameras[cam1], newCameras[cam2], 0);
	e1 = eg.computeClosest3DPoint(imageCoordinates1[1]/2.7, imageCoordinates2[1]/2.7, newCameras[cam1], newCameras[cam2], 0);
	
	//line2
	s2 = eg.computeClosest3DPoint(imageCoordinates1[2]/2.7, imageCoordinates2[2]/2.7, newCameras[cam1], newCameras[cam2], 0);
	e2 = eg.computeClosest3DPoint(imageCoordinates1[3]/2.7, imageCoordinates2[3]/2.7, newCameras[cam1], newCameras[cam2], 0);
		
	printf("newScale:%f\n", (s1cur-e1cur).length() / (s1 - e1).length());
		
	//rotation-translation
	//cam1 has rotation R1 and center1, a point x in camera coordinates maps to Y in world coordinates
	//we want the camera to have a rotation R2 and center C2 so that the same point in camera coordinates,
	//x, maps to world coordinate Z 
	Matrix4f extMat = expand3To4(KInvMats[cam1]) * newCameras[cam1].getProjectionMatrix(0);
	Vec3f center = shrink4To3(newCameras[cam1].getCenter());
		
	vector<Vec3f> points;
	points.push_back(shrink4To3(extMat*expand3To4(s1)));
	points.push_back(shrink4To3(extMat*expand3To4(e1)));
	points.push_back(shrink4To3(extMat*expand3To4(s2)));
	points.push_back(shrink4To3(extMat*expand3To4(e2)));
		
	Matrix3f rot;
	Vec3f translation;
	GeometryUtils::findRigidTransformation(points, refPoints, rot, translation);
		
	printf("rotation:");
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			printf("%f ", rot[j][i]);
		}
		printf("\n");
	}
	
	printf("translation:%f %f %f\n", translation[0], translation[1], translation[2]);
		
	printf("residuals:\n");
	for(int i=0; i<points.size(); i++)
	{
		printf("%f\n", (refPoints[i] - (rot*points[i]+translation)).length());
	}
		
	//pc->rotateData(rot*shrink4To3(extMat));
	//pc->translateData(translation+rot*shrink4To3(extMat[3]));
	pc->writeMesh("/Users/ceylan/Desktop/deneme.ply");
		
	//transform the reference point so that imageIndex1 has the new extrinsic parameters
	rot = rot.transpose();
	Matrix4f transformation = expand3To4((shrink4To3(extMat)).transpose()*rot);	
	Vec3f t = center - shrink4To3(transformation) * translation;
	transformation[3][0] = t[0]; transformation[3][1] = t[1]; transformation[3][2] = t[2];
		
	string calibResultPath = mainDirectory + "/comparisonJiang2012/newCam";
	mkdir(calibResultPath.c_str(), 0770);
	
	for(int i=0;i<newCameras.size();i++)
	{
        //scale
		if(i!=cam1 && i!=cam2)
		{
			newCenter = shrink4To3(newCameras[i].getCenter()) * scale;
			newCenter = -shrink4To3(newCameras[i].getProjectionMatrix(0)) * newCenter;
			scaleMat = newCameras[i].getProjectionMatrix(0);
			scaleMat[3][0] = newCenter[0]; scaleMat[3][1] = newCenter[1]; scaleMat[3][2] = newCenter[2];
			newCameras[i].setProjectionMatrix(scaleMat);
		}
			
		Matrix4f extMat = expand3To4(KInvMats[i]) * newCameras[i].getProjectionMatrix(0);
		extMat = extMat * transformation;
		//newCameras[i].setProjectionMatrix(expand3To4(KMats[i])*extMat);
		
		Matrix4f P = newCameras[i].getProjectionMatrix(0);
		char buf[1024];
		sprintf(buf, "%s/%08d.txt", calibResultPath.c_str(), i);
		
		ofstream foutCam(buf);
		foutCam << "CONTOUR" << endl;
		foutCam << P[0][0] << " " << P[1][0] << " " << P[2][0] << " " << P[3][0] << endl;
		foutCam << P[0][1] << " " << P[1][1] << " " << P[2][1] << " " << P[3][1] << endl;
		foutCam << P[0][2] << " " << P[1][2] << " " << P[2][2] << " " << P[3][2] << endl;
	}
	
	Mesh3D camMesh;
	
	for(int i=0; i<newCameras.size(); i++)
	{
		Vec3f center = shrink4To3(newCameras[i].getCenter());
		Vec3f ll = newCameras[i].unprojectPixelwithDepth(Vec3f(0.0, 800/*imageSizes[i][1]*/, 1), 0);
		Vec3f lr = newCameras[i].unprojectPixelwithDepth(Vec3f(1067/*imageSizes[i][0]*/, 800/*imageSizes[i][1]*/, 1), 0);
		Vec3f ul = newCameras[i].unprojectPixelwithDepth(Vec3f(0.0, 0.0, 1), 0);
		Vec3f ur = newCameras[i].unprojectPixelwithDepth(Vec3f(1067/*imageSizes[i][0]*/, 0.0, 1), 0);
		
		Vec3f dir;
		dir = ll - center;
		dir.normalize();
		ll = center + dir*1.5;
		
		dir = lr - center;
		dir.normalize();
		lr = center + dir*1.5;
		
		dir = ul - center;
		dir.normalize();
		ul = center + dir*1.5;
		
		dir = ur - center;
		dir.normalize();
		ur = center + dir*1.5;
		
		camMesh.addPoint(center);
		camMesh.addPoint(ll);
		camMesh.addPoint(lr);
		camMesh.addPoint(ul);
		camMesh.addPoint(ur);
		
		camMesh.addFace(i*5, i*5+3, i*5+1);
		camMesh.addFace(i*5, i*5+1, i*5+2);
		camMesh.addFace(i*5, i*5+2, i*5+4);
		camMesh.addFace(i*5, i*5+4, i*5+3);
		
		camMesh.addFace(i*5+1, i*5+3, i*5+4);
		camMesh.addFace(i*5+4, i*5+2, i*5+1);
	}
	string filename = "/Users/ceylan/Desktop/cameras.obj";
	camMesh.writeMesh(filename.c_str());
	
}

void BundlerManager::read3DGrids()
{
	string filename = mainDirectory + "/matches/grid.txt";
	ifstream gridFile(filename.c_str(), ios::in);
	
	if(!gridFile)
	{
		printf("3D grid file not found!\n");
		return;
	}
	
	string repFilename = mainDirectory + "/matches/repetitionCells.txt";
	ifstream repFile(repFilename.c_str(), ios::in);
	
	if(!repFile)
	{
		printf("Repetition file not found!\n");
		return;
	}
	
	if(sparsePC == NULL)
	{
		printf("No sparse point cloud found!\n");
		return;
	}
	
	EpipolarGeometry eg;
	
	Vec3f lowerPt, upperPt;
	sparsePC->getBoundingBox(lowerPt, upperPt);
	float diag = (upperPt-lowerPt).length();
	float threshold = diag / 100;
	
	int noImages, noGrids;
	repFile >> noImages >> noGrids;
	
	gridFile >> noGrids;
	
	grids3D.resize(noGrids);
	
	for(int g=0; g<noGrids; g++)
	{
		int tmpImg, tmpPlane;
		Vec2f tmpCenter;
		Vec2i tmpSize;
		
		int templateId;
		
		gridFile >> templateId;
		
		int noCols, noRows;
		Vec3f pt, horTrans, verTrans;
		
		repFile >> noRows >> noCols;
		repFile >> tmpImg >> tmpPlane >> tmpSize[0] >> tmpSize[1] >> tmpCenter[0] >> tmpCenter[1];
		
		gridFile >> pt[0] >> pt[1] >> pt[2];
		gridFile >> horTrans[0] >> horTrans[1] >> horTrans[2];
		gridFile >> verTrans[0] >> verTrans[1] >> verTrans[2];
		grids3D[g] = TranslationalGrid3D(noCols, noRows);
		grids3D[g].setTransformations(horTrans, verTrans);
		grids3D[g].setTemplateInfo(tmpImg, tmpPlane, tmpCenter, tmpSize);
		grids3D[g].fillElementCenters(pt);
		
		vector<vector<Vec2f> > repMask;
		vector<int> planeIndices;
		vector<int> gridIndices;
		
		repMask.resize(noImages);
		planeIndices.resize(noImages, -1);
		gridIndices.resize(noImages, -1);
		
		for(int i=0; i<noImages; i++)
		{
			printf("img %d\n", i);
			for(int c=0; c<noCols; c++)
			{
				for(int r=0; r<noRows; r++)
				{
					int plane, grid;
					Vec2f pos;
					repFile >> plane >> grid >> pos[0] >> pos[1];
					printf("g:%d %d %f %f\n", plane, grid, pos[0], pos[1]);
					repMask[i].push_back(pos);
					if(plane != -1)
						planeIndices[i] = plane;
					if(grid != -1)
						gridIndices[i] = grid;
					if(plane != -1 && grid != -1)
						grids3D[g].addVisibleImage(i);
				}
			}
		}
		
		if(noRows == 1)
		{
			int index1, index2;
				
			bool found = false;
			int colIndex;
			Vec2f refPixel, neigPixel;
			for(int img1=0; img1<noImages; img1++)
			{
				for(int img2=0; img2<noImages; img2++)
				{
					if(img2==img1)
						continue;
						
					for(int c=0; c<noCols; c++)
					{
						for(int r=0; r<noRows; r++)
						{
							if(repMask[img1][c*noRows+r] != Vec2f(0.0, 0.0) && repMask[img2][c*noRows+r] != Vec2f(0.0, 0.0))
							{
								colIndex = c;
								found = true;
								index1 = img1;
								index2 = img2;
								refPixel = repMask[img1][c*noRows+r];
								neigPixel = repMask[img2][c*noRows+r];
								break;
							}
						}
						if(found)
							break;
					}
					if(found)
						break;
				}
				if (found)
					break;
			}
			
			refPixel -= Vec2f(0.0, repeatingGrids[templateId][index1][planeIndices[index1]][gridIndices[index1]].gridCells[0].size[1]/2.0);
			
			neigPixel -= Vec2f(0.0, repeatingGrids[templateId][index2][planeIndices[index2]][gridIndices[index2]].gridCells[0].size[1]/2.0);
			
			Vec3f refPixel2 = rectifyingHomographies[index1][planeIndices[index1]] * expand2To3(refPixel);
			refPixel2[0] = refPixel2[0]/refPixel2[2];
			refPixel2[1] = refPixel2[1]/refPixel2[2];
				
			Vec3f neigPixel2 = rectifyingHomographies[index2][planeIndices[index2]] * expand2To3(neigPixel);
			neigPixel2[0] = neigPixel2[0]/neigPixel2[2];
			neigPixel2[1] = neigPixel2[1]/neigPixel2[2];
				
			Vec3f p1 = eg.computeClosest3DPoint(shrink3To2(refPixel2), shrink3To2(neigPixel2), cameras[index1], cameras[index2], 0);
				
				
			Vec3f p2 = pt+horTrans*colIndex;
			
			verTrans = p2 - p1;
		}
		
		Vec3f h = horTrans;
		h.normalize();
		Vec3f v = verTrans;
		v.normalize();
		Vec3f n = cross(h, v);
		n.normalize();
		float d = 0.0 - dot(n, pt);
		printf("n:%f %f %f, d:%f\n", n[0], n[1], n[2], d);
		Plane pln = Plane(n,d);
		pln.setAxis(h, v);
		pln.setPlaneCentroid(pt);
		float minX, minY, maxX, maxY;
		
		minX = dot(h, pt-pln.getPlaneCentroid());
		minY = dot(v, pt-pln.getPlaneCentroid());
		
		maxX = dot(h, pt+horTrans*noCols+verTrans*noRows-pln.getPlaneCentroid());
		maxY = dot(v, pt+horTrans*noCols+verTrans*noRows-pln.getPlaneCentroid());
		
		vector<Vec3f> inliers;
		for(int p=0; p<sparsePC->getNoPoints(); p++)
		{
			float dist = dot(n, sparsePC->getVertexPosition(p)) + d;
			//printf("dist:%f\n", dist);
			
			if(abs(dist) < threshold)
			{
				inliers.push_back(sparsePC->getVertexPosition(p));
				float xD = dot(h, sparsePC->getVertexPosition(p)-pt);
				float yD = dot(v, sparsePC->getVertexPosition(p)-pt);
				if(xD < minX)
					minX = xD;
				else if(xD > maxX)
					maxX = xD;
				if(yD < minY)
					minY = yD;
				else if(yD > maxY)
					maxY = yD;
			}
		}
		
		planeInliers.push_back(inliers);
		
		PointCloud *tmp = new PointCloud();
		float xIncr = (maxX-minX)/1000.0;
		float yIncr = (maxY-minY)/1000.0;
		for(int x=0; x<1000; x++)
		{
			for(int y=0; y<1000; y++)
			{
				Vec3f pos = pln.getPlaneCentroid()+h*(minX+xIncr*x)+v*(minY+yIncr*y);
				tmp->addPoint(pos);
			}
		}
		
		tmp->writeMesh("/Users/ceylan/Desktop/pl1.ply");
		
		pln.setBoundingBox(pln.getPlaneCentroid()+h*minX+v*minY, pln.getPlaneCentroid()+h*maxX+v*maxY);
		grids3D[g].setGridPlane(pln);
		
		bool samePlane = false;
		for(int p=0; p<planes.size(); p++)
		{
			if(acos(dot(n, planes[p].getPlaneNormal()))*180.0 / M_PI < 10.0)
			{
				if(abs(d-planes[p].getPlaneDistance()) < threshold)
				{
					samePlane = true;
					break;
				}
			}
		}
		
		if(!samePlane)
			planes.push_back(pln);
	}
	
}

void BundlerManager::fitPlanesTo3DData()
{
	if(sparsePC == NULL)
		return ;
	int noPoints = sparsePC->getNoPoints();
	Vec3f lowerPt, upperPt;
	sparsePC->getBoundingBox(lowerPt, upperPt);
	float diag = (upperPt-lowerPt).length();
	float threshold = diag / 100.0;
	
	vector<bool> selectedMask;
	selectedMask.resize(noPoints, false);
	
	//ransac loop
	int noPtsUsed = 3;
	int numTrials = 512;
	
	bool stop = false;
	
	string filename = mainDirectory + "/matches/grid.txt";
	ifstream gridFile(filename.c_str(), ios::in);
	
	if(!gridFile)
	{
		printf("3D grid file not found!\n");
		return;
	}
	
	string repFilename = mainDirectory + "/matches/repetitionCells.txt";
	ifstream repFile(repFilename.c_str(), ios::in);
	
	if(!repFile)
	{
		printf("Repetition file not found!\n");
		return;
	}
	
	int noImages, noGrids;
	repFile >> noImages >> noGrids;
	
	gridFile >> noGrids;
	
	grids3D.resize(noGrids);
	
	for(int g=0; g<noGrids; g++)
	{
		int tmpImg, tmpPlane;
		Vec2f tmpCenter;
		Vec2i tmpSize;
	
		int templateId;
	
		gridFile >> templateId;
	
		int noCols, noRows;
		Vec3f pt, horTrans, verTrans;
	
		repFile >> noRows >> noCols;
		repFile >> tmpImg >> tmpPlane >> tmpSize[0] >> tmpSize[1] >> tmpCenter[0] >> tmpCenter[1];
	
		gridFile >> pt[0] >> pt[1] >> pt[2];
		gridFile >> horTrans[0] >> horTrans[1] >> horTrans[2];
		gridFile >> verTrans[0] >> verTrans[1] >> verTrans[2];
		grids3D[g] = TranslationalGrid3D(noCols, noRows);
		grids3D[g].setTransformations(horTrans, verTrans);
		grids3D[g].setTemplateInfo(tmpImg, tmpPlane, tmpCenter, tmpSize);
		grids3D[g].fillElementCenters(pt);
		
		vector<vector<Vec2f> > repMask;
		vector<int> planeIndices;
		vector<int> gridIndices;
		
		repMask.resize(noImages);
		planeIndices.resize(noImages, -1);
		gridIndices.resize(noImages, -1);
		
		for(int i=0; i<noImages; i++)
		{
			printf("img %d\n", i);
			for(int c=0; c<noCols; c++)
			{
				for(int r=0; r<noRows; r++)
				{
					int plane, grid;
					Vec2f pos;
					repFile >> plane >> grid >> pos[0] >> pos[1];
					printf("g:%d %d %f %f\n", plane, grid, pos[0], pos[1]);
					repMask[i].push_back(pos);
					if(plane != -1)
						planeIndices[i] = plane;
					if(grid != -1)
						gridIndices[i] = grid;
					if(plane != -1 && grid != -1)
						grids3D[g].addVisibleImage(i);
				}
			}
		}
	}
	
		
	while(!stop)
	{
		stop = true;
		
		vector<Vec3f> normals;
		vector<float> distances;
		vector<int> inliers;
		int maxInlierIndex = -1;
		
		for(int i=0; i<numTrials; i++)
		{
			int *idxs = new int[noPtsUsed];
			int round = 0;
		
			// Sample 3 random correspondences
			int j;
			for (j = 0; j < noPtsUsed; j++) 
			{
				int reselect = 0;
			
				if (round == 1000)
					break;
			
				int idx = rand() % noPoints;
			
				if(selectedMask[idx] == true)
				{
					reselect = 1;
				}
				
				// Make sure we didn't sample this index yet
				for (int k = 0; k < j; k++) 
				{
					if (idx == idxs[k])
					{
						reselect = 1;
						break;
					}
				}
			
				if (reselect) 
				{
					round++;
					j--;
					continue;
				}
			
				idxs[j] = idx;
			}
		
			if(j != noPtsUsed)
				break;
			
			//find plane parameters
			Vec3f v1 = sparsePC->getVertexPosition(idxs[2]) - sparsePC->getVertexPosition(idxs[0]);
			Vec3f v2 = sparsePC->getVertexPosition(idxs[1]) - sparsePC->getVertexPosition(idxs[0]);
			Vec3f n = cross(v1.normalize(), v2.normalize());
			n.normalize();
			if(n == Vec3f(0.0, 0.0, 0.0))
				continue;
			
			float d = 0.0 - dot(n, sparsePC->getVertexPosition(idxs[0]));
		
			normals.push_back(n);
			distances.push_back(d);
		
			//find inliers
			int noInliers = 0;
			for(int p=0; p<noPoints; p++)
			{
				if(selectedMask[p])
					continue;
				float dist = dot(n, sparsePC->getVertexPosition(p)) + d;
				if(abs(dist) < threshold)
					noInliers += 1;
			}
			inliers.push_back(noInliers);
			if(maxInlierIndex == -1)
				maxInlierIndex = inliers.size()-1;
			else if(noInliers > inliers[maxInlierIndex])
				maxInlierIndex = inliers.size()-1;
		}
		printf("maxInliers:%d\n", inliers[maxInlierIndex]);
		if(inliers[maxInlierIndex] > 1000)
		{
			printf("*******plane******\n");
			int count = 0;
			
			vector<Vec3f> inlierPts;
			
			for(int p=0; p<noPoints; p++)
			{
				if(selectedMask[p])
					continue;
				float dist = dot(normals[maxInlierIndex], sparsePC->getVertexPosition(p)) + distances[maxInlierIndex];
				if(abs(dist) < threshold)
				{
					selectedMask[p] = true;
					//project to plane
					inlierPts.push_back(sparsePC->getVertexPosition(p));
					count += 1;
				}
			}
			
			Vec3f n, centroid;
			GeometryUtils::projectDataToPlane(inlierPts, n, centroid);
			float d = 0.0 - dot(n,centroid);

			Vec3f xAxis, yAxis;
			for(int g=0; g<noGrids; g++)
			{
				Vec3f pt = grids3D[g].getElementCenter(0,0);
				if(abs(dot(n, pt) + d) < threshold)
				{
					printf("grid:%d\n", g);
					//find x and y axes
					Plane pln(n,d);
					planes.push_back(pln);
					int templateImg, templatePlane;
					Vec2f templateCenter;
					Vec2i templateSize;
					grids3D[g].getTemplateInfo(templateImg,templatePlane,templateCenter,templateSize);
					//x-axis
					Vec3f p1 = rectifyingHomographies[templateImg][templatePlane]*Vec3f(templateCenter[0]-templateSize[0]/2.0, templateCenter[1],1.0);
					Vec3f p2 = rectifyingHomographies[templateImg][templatePlane]*Vec3f(templateCenter[0]+templateSize[0]/2.0, templateCenter[1],1.0);
					p1[0] /= p1[2]; p1[1] /= p1[2]; p1[2] /= p1[2];
					p2[0] /= p2[2]; p2[1] /= p2[2]; p2[2] /= p2[2];
					Vec3f rayDir1 = cameras[templateImg].unprojectPixelwithDepth(p1, 0) - shrink4To3(cameras[templateImg].getCenter());
					Vec3f rayDir2 = cameras[templateImg].unprojectPixelwithDepth(p2, 0) - shrink4To3(cameras[templateImg].getCenter());
					Vec3f projP1, projP2;
					pln.findPlaneRayIntersection(shrink4To3(cameras[templateImg].getCenter()), rayDir1, projP1);
					pln.findPlaneRayIntersection(shrink4To3(cameras[templateImg].getCenter()), rayDir2, projP2);
					xAxis = (projP2-projP1);
					xAxis.normalize();
					//y-axis
					p1 = rectifyingHomographies[templateImg][templatePlane]*Vec3f(templateCenter[0], templateCenter[1]-templateSize[1]/2.0,1.0);
					p2 = rectifyingHomographies[templateImg][templatePlane]*Vec3f(templateCenter[0], templateCenter[1]+templateSize[1]/2.0,1.0);
					p1[0] /= p1[2]; p1[1] /= p1[2]; p1[2] /= p1[2];
					p2[0] /= p2[2]; p2[1] /= p2[2]; p2[2] /= p2[2];
					rayDir1 = cameras[templateImg].unprojectPixelwithDepth(p1, 0) - shrink4To3(cameras[templateImg].getCenter());
					rayDir2 = cameras[templateImg].unprojectPixelwithDepth(p2, 0) - shrink4To3(cameras[templateImg].getCenter());
					pln.findPlaneRayIntersection(shrink4To3(cameras[templateImg].getCenter()), rayDir1, projP1);
					pln.findPlaneRayIntersection(shrink4To3(cameras[templateImg].getCenter()), rayDir2, projP2);
					yAxis = projP2-projP1;
					yAxis.normalize();
					
					xAxis = cross(yAxis, n).normalize();
					yAxis = cross(n, xAxis).normalize();
					
					planes[planes.size()-1].setPlaneCentroid(centroid);
					planes[planes.size()-1].setAxis(xAxis, yAxis);
			
					planeInliers.push_back(inlierPts);
			
					printf("plane with n:%f, %f, %f dist:%f\n",normals[maxInlierIndex][0], normals[maxInlierIndex][1], normals[maxInlierIndex][2], distances[maxInlierIndex]);
					printf("%d points added\n", count);
					stop = false;
			
					float minX, minY, maxX, maxY;
		
					PointCloud *tmp = new PointCloud();
					for(int p=0; p<sparsePC->getNoPoints(); p++)
					{
						float dist = dot(n, sparsePC->getVertexPosition(p)) + d;
						if(abs(dist) < threshold)
						{
							tmp->addPoint(sparsePC->getVertexPosition(p));
							float xD = dot(xAxis, sparsePC->getVertexPosition(p)-centroid);
							float yD = dot(yAxis, sparsePC->getVertexPosition(p)-centroid);
							if(p==0)
							{
								minX = xD; maxX = xD;
								minY = yD; maxY = yD;
							}
							else
							{
								if(xD < minX)
									minX = xD;
								else if(xD > maxX)
									maxX = xD;
								if(yD < minY)
									minY = yD;
								else if(yD > maxY)
									maxY = yD;
							}
						}
					}
					
					float xIncr = (maxX-minX)/1000.0;
					float yIncr = (maxY-minY)/1000.0;
					for(int x=0; x<1000; x++)
					{
						for(int y=0; y<1000; y++)
						{
							tmp->addPoint(centroid+xAxis*(minX+xIncr*x)+yAxis*(minY+yIncr*y));
						}
					}
					tmp->writeMesh("/Users/ceylan/Desktop/pl1.ply");
					planes[planes.size()-1].setBoundingBox(planes[planes.size()-1].getPlaneCentroid()+xAxis*minX+yAxis*minY, planes[planes.size()-1].getPlaneCentroid()+xAxis*maxX+yAxis*maxY);
					grids3D[g].setGridPlane(planes[planes.size()-1]);
					
					int rows = grids3D[g].getNoRows();
					int cols = grids3D[g].getNoColumns();
					for(int r=0; r<rows; r++)
					{
						for(int c=0; c<cols; c++)
						{
							Vec3f pt = grids3D[g].getElementCenter(r,c);
							pt = pt - n*dot(n,pt-centroid);
							grids3D[g].setElementCenter(r,c, pt);
						}
					}
				}
			}
		}
	}
	
	/*for(int i=0; i<planes.size(); i++)
	{
		vector<Vec3f> directions;
		for(int j=0; j<planes.size(); j++)
		{
			if(i==j)
				continue;
			
			float dotProd = dot(planes[i].getPlaneNormal(), planes[j].getPlaneNormal());
			if(abs(dotProd) < 0.5)
			{
				Vec3f dir = cross(planes[i].getPlaneNormal(), planes[j].getPlaneNormal());
				dir.normalize();
				directions.push_back(dir);
				printf("dir:%f %f %f\n", dir[0], dir[1], dir[2]);
			}
		}
		float minAngleDiff;
		int repId = -1;
		for(int d=0; d<directions.size(); d++)
		{
			float diff = 0.0;
			for(int e=0; e<directions.size(); e++)
			{
				if(d==e)
					continue;
				diff += abs(dot(directions[d], directions[e]));
			}
			if(repId == -1)
			{
				repId = d;
				minAngleDiff = diff;
			}
			else if(minAngleDiff < diff)
			{
				repId = d;
				minAngleDiff = diff;
			}
		}
		
		Vec3f yAxis = directions[repId];
		Vec3f xAxis = cross(planes[i].getPlaneNormal(), yAxis);
		planes[i].setAxis(xAxis, yAxis);
		
		float minX, minY, maxX, maxY;
		for(int p=0; p<planeInliers[i].size(); p++)
		{
			float xD = dot(xAxis, planeInliers[i][p]-planes[i].getPlaneCentroid());
			float yD = dot(yAxis, planeInliers[i][p]-planes[i].getPlaneCentroid());
			if(p==0)
			{
				minX = xD; maxX = xD;
				minY = yD; maxY = yD;
			}
			else
			{
				if(xD < minX)
					minX = xD;
				else if(xD > maxX)
					maxX = xD;
				if(yD < minY)
					minY = yD;
				else if(yD > maxY)
					maxY = yD;
			}
		}
		
		PointCloud *tmp = new PointCloud();
		
		tmp->addPoint(planes[i].getPlaneCentroid()+xAxis);
		tmp->addPoint(planes[i].getPlaneCentroid()+yAxis);
		for(int x=(int)minX; x<(int)maxX; x++)
		{
			for(int y=(int)minY; y<(int)maxY; y++)
			{
				Vec3f pos = planes[i].getPlaneCentroid()+xAxis*x+yAxis*y;
				tmp->addPoint(pos);
			}
		}
		
		tmp->writeMesh("/Users/ceylan/Desktop/pl1.ply");
		planes[i].setBoundingBox(planes[i].getPlaneCentroid()+xAxis*minX+yAxis*minY, planes[i].getPlaneCentroid()+xAxis*maxX+yAxis*maxY);
	}*/
	
}

void BundlerManager::rectifyWith3DPlane()
{		
	int noImages = imageFileNames.size();
	
	for(int index=30; index<noImages; index++)
	{
		//find the planes
		vector<int> inlierCount;
		inlierCount.resize(planes.size(), 0);
	
		for(int i=0; i<grids3D.size(); i++)
		{
			for(int p=0; p<planeInliers[i].size(); p++)
			{
				//project to the image
				Vec3f proj = cameras[index].project(expand3To4(planeInliers[i][p]), 0);
				if(proj[0] >= 0.0 && proj[0] < imageSizes[index][0] && proj[1] >= 0.0 && proj[1] < imageSizes[index][1])
				{
					inlierCount[i] = inlierCount[i] + 1;
				}
			}
		}
	
		Img currentImg;
		currentImg.read(imageFileNames[index]);
	
		for(int i=0; i<grids3D.size(); i++)
		{
			char imgFile[1024];
			sprintf(imgFile, "%08d_%d_final.%s", index, i, "jpg");
			
			string rectFileName = mainDirectory + "/rectification/" + imgFile;
			ifstream rectFile(rectFileName.c_str(), ios::in);
			if(rectFile)
				continue;
			
			int width = rectImgWidth;
			int height = rectImgHeight;
			
			if(inlierCount[i] > 50)
			{
				Img rectifiedImage;
				rectifiedImage.resize(width, height, Color(0,0,0));
				Vec3f minCoord, maxCoord;
				Plane pln = grids3D[i].getGridPlane();
				pln.getBoundingBox(minCoord, maxCoord);
				Vec3f xAxis = pln.getPlaneXAxis();
				Vec3f yAxis = pln.getPlaneYAxis();
				float xIncr = dot(xAxis, (maxCoord-minCoord)) / width;
				float yIncr =  dot(yAxis, (maxCoord-minCoord)) / height;
				PointCloud *tmp = new PointCloud();
				
				for(int x=0; x<width; x++)
				{
					for(int y=0; y<height; y++)
					{
						Vec3f pt = minCoord + xAxis*xIncr*x + yAxis*yIncr*y;
						tmp->addPoint(pt);
						Vec3f proj = cameras[index].project(expand3To4(pt), 0);
						if(proj[0] >= 0.0 && proj[0] < imageSizes[index][0] && proj[1] >= 0.0 && proj[1] < imageSizes[index][1])
						{
							rectifiedImage.setColor(x, y, currentImg(proj[0], proj[1]));
						}
					}
				}
				rectifiedImage.write(rectFileName);
				tmp->writeMesh("/Users/ceylan/Desktop/p1.ply");
			}
			printf("plane %d:%d\n", i, inlierCount[i]);
		}
	}
}

void BundlerManager::updateGridCells()
{
	int noImages = imageFileNames.size();
	rectifyingHomographies.resize(noImages);
	noRectifiedImages.resize(noImages, 0);
	
	RepetitionFinder rf;
	float threshold = 0.6;
	
	for(int i=0; i<noImages; i++)
	{
		rectifyingHomographies[i].resize(2);
		
		char homFileName[1024];
		sprintf(homFileName, "homography%d_0.txt", i);
		string homographyFileName = mainDirectory + "/rectification/" + homFileName;
		
		char newHomFileName[1024];
		sprintf(newHomFileName, "homography%d_1.txt", i);
		string newHomographyFileName = mainDirectory + "/rectification/" + newHomFileName;
		
		ifstream homFile(homographyFileName.c_str());
		if(homFile)
		{
			Matrix3f H;
			homFile >> H[0][0] >> H[1][0] >> H[2][0];
			homFile >> H[0][1] >> H[1][1] >> H[2][1];
			homFile >> H[0][2] >> H[1][2] >> H[2][2];
			homFile.close();
			rectifyingHomographies[i][0] = H;
			noRectifiedImages[i] = noRectifiedImages[i]+1;
		}
			
		ifstream newHomFile(newHomographyFileName.c_str(), ios::in);
		if(newHomFile)
		{
			Matrix3f HSecond;
				
			newHomFile >> HSecond[0][0] >> HSecond[1][0] >> HSecond[2][0];
			newHomFile >> HSecond[0][1] >> HSecond[1][1] >> HSecond[2][1];
			newHomFile >> HSecond[0][2] >> HSecond[1][2] >> HSecond[2][2];
			newHomFile.close();
			rectifyingHomographies[i][1] = HSecond;
			noRectifiedImages[i] = noRectifiedImages[i]+1;
		}
	}
	
	int noGrids = grids3D.size();
	
	int count = 0;
	for(int g=0; g<noGrids; g++)
	{
		
		//find the template cell
		int tmpImg, tmpPlane;
		Vec2f tmpCenter;
		Vec2i tmpSize;
		
		grids3D[g].getTemplateInfo(tmpImg, tmpPlane, tmpCenter, tmpSize);
		int w = tmpSize[0]; int h = tmpSize[1];
		vector<Vec3f> templateCoordinates;
		Vec3f raydir;
		
		Vec3f p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0] - w/2, tmpCenter[1] - h/2, 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		templateCoordinates.push_back(p);
		
		p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0] + w/2, tmpCenter[1] - h/2, 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		templateCoordinates.push_back(p);
		
		p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0] + w/2, tmpCenter[1] + h/2, 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		templateCoordinates.push_back(p);
		
		p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0] - w/2, tmpCenter[1] + h/2, 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		templateCoordinates.push_back(p);
		
		Vec3f minCoord, maxCoord;
		Vec3f xAxis, yAxis;
		Plane pln = grids3D[g].getGridPlane();
		pln.getBoundingBox(minCoord, maxCoord);
		xAxis = pln.getPlaneXAxis(); yAxis = pln.getPlaneYAxis();
		
		float xIncr = dot(xAxis, (maxCoord-minCoord)) / rectImgWidth;
		float yIncr =  dot(yAxis, (maxCoord-minCoord)) / rectImgHeight;
		
		float minX, maxX, minY, maxY;
		for(int i=0; i<templateCoordinates.size(); i++)
		{
			float xD = dot(xAxis, templateCoordinates[i] - minCoord) / xIncr;
			float yD = dot(yAxis, templateCoordinates[i] - minCoord) / yIncr;
			
			if(i==0)
			{
				minX = xD; maxX = xD; minY = yD; maxY = yD;
			}
			else 
			{
				if(xD < minX)
					minX = xD;
				else if(xD>maxX)
					maxX = xD;
				if(yD < minY)
					minY = yD;
				else if(yD > maxY)
					maxY = yD;
			}
		}
		
		minX = floor(minX); minY = floor(minY);
		maxX = floor(maxX); maxY = floor(maxY);
		
		Img templateImg;
		templateImg.resize(maxX-minX, maxY-minY);
		char imgFile[1024];
		sprintf(imgFile, "%08d_%d_final.%s", tmpImg, g, "jpg");
		
		string rectFileName = mainDirectory + "/rectification/" + imgFile;
		Img rectImg;
		rectImg.read(rectFileName);
		
		for(int x=minX; x<maxX; x++)
		{
			for(int y=minY; y<maxY; y++)
			{
				templateImg.setColor(x-minX, y-minY, rectImg(x, y));
			}
		}
		
		templateImg.write("/Users/ceylan/Desktop/tmp.png");
		
		Vec2i searchSize((maxX-minX)*2, (maxY-minY)*2);
		
		for(int i=0; i<noImages; i++)
		{
			printf("*****img %d*****\n", i);
			char imgFile[1024];
			sprintf(imgFile, "%08d_%d_final.%s", i, g, "jpg");
			
			string rectFileName = mainDirectory + "/rectification/" + imgFile;
			Img rectImg;
			rectImg.read(rectFileName);
			
			for(int r=0; r<grids3D[g].getNoRows(); r++)
			{
				for(int c=0; c<grids3D[g].getNoColumns(); c++)
				{
					Vec3f center = grids3D[g].getElementCenter(r, c);
					
					Vec3f pos = cameras[i].project(expand3To4(center), 0);
					if(pos[0] >= 0.0 && pos[0] < imageSizes[i][0] && pos[1] >= 0.0 && pos[1] < imageSizes[i][1])
					{
						printf("pos:%f %f\n", pos[0], pos[1]);
						keyInfo[i][count + r*grids3D[g].getNoColumns() + c].x = pos[0];
						keyInfo[i][count + r*grids3D[g].getNoColumns() + c].y = pos[1];
					}
					
					/*float xD = dot(xAxis, center - minCoord) / xIncr;
					float yD = dot(yAxis, center - minCoord) / yIncr;
					
					
					Vec2i offset(xD-searchSize[0]/2, yD-searchSize[1]/2);
					vector<Vec2f> tempMatches;
					vector<float> scores;
			
					rf.templateMatchingWithFFTCorrelation(&templateImg, &rectImg, offset, searchSize, tempMatches, scores, threshold);
					if(tempMatches.size() > 0)
					{
						Vec3f pos = minCoord + xAxis*xIncr*tempMatches[0][0] + yAxis*yIncr*tempMatches[0][1];
						pos = cameras[i].project(expand3To4(pos), 0);
						if(pos[0] >= 0.0 && pos[0] < imageSizes[i][0] && pos[1] >= 0.0 && pos[1] < imageSizes[i][1])
						{
							printf("pos:%f %f\n", pos[0], pos[1]);
							keyInfo[i][count + r*grids3D[g].getNoColumns() + c].x = pos[0];
							keyInfo[i][count + r*grids3D[g].getNoColumns() + c].x = pos[1];
						}
					}*/
				}
			}
		}
		
		count += grids3D[g].getNoColumns() * grids3D[g].getNoRows();
	}
	
	for(int i=0; i<noImages; i++)
	{
		char keyFile[1024];
		sprintf(keyFile, "%08d.%s", i, "key");
		string filename = mainDirectory + "/matches/finalMatches/" + keyFile;
		ofstream fout(filename.c_str());
		fout << keyInfo[i].size() << " 128\n";
		
		for(int r=0; r<keyInfo[i].size(); r++)
		{
			if(r < count)
				fout << keyInfo[i][r].y << " " << keyInfo[i][r].x << " 0.0 0.0 1\n";
			else
				fout << keyInfo[i][r].y << " " << keyInfo[i][r].x << " 0.0 0.0 0\n";
			for(int d=0; d<7; d++)
			{
				if(d==6)
					fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
				else
					fout << "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n";
			}
		}
		
		fout.close();
	}
}

bool BundlerManager::form3DGrids()
{
	int noGrids = grids3D.size();
	
	int noImages = imageFileNames.size();
	
	computeEdges();
	
	for(int i=0; i<noImages; i++)
	{
		rectifyingHomographies[i].resize(2);
		
		char homFileName[1024];
		sprintf(homFileName, "homography%d_0.txt", i);
		string homographyFileName = mainDirectory + "/rectification/" + homFileName;
		
		char newHomFileName[1024];
		sprintf(newHomFileName, "homography%d_1.txt", i);
		string newHomographyFileName = mainDirectory + "/rectification/" + newHomFileName;
		
		ifstream homFile(homographyFileName.c_str());
		if(homFile)
		{
			Matrix3f H;
			homFile >> H[0][0] >> H[1][0] >> H[2][0];
			homFile >> H[0][1] >> H[1][1] >> H[2][1];
			homFile >> H[0][2] >> H[1][2] >> H[2][2];
			homFile.close();
			rectifyingHomographies[i][0] = H;
			noRectifiedImages[i] = noRectifiedImages[i]+1;
		}
		
		ifstream newHomFile(newHomographyFileName.c_str(), ios::in);
		if(newHomFile)
		{
			Matrix3f HSecond;
			
			newHomFile >> HSecond[0][0] >> HSecond[1][0] >> HSecond[2][0];
			newHomFile >> HSecond[0][1] >> HSecond[1][1] >> HSecond[2][1];
			newHomFile >> HSecond[0][2] >> HSecond[1][2] >> HSecond[2][2];
			newHomFile.close();
			rectifyingHomographies[i][1] = HSecond;
			noRectifiedImages[i] = noRectifiedImages[i]+1;
		}
	}
	
	float score;
	int boundaryThresh = 8;
	int index;
	
	for(int g=0; g<noGrids; g++)
	{
		int tmpImg, tmpPlane;
		Vec2f tmpCenter;
		Vec2i tmpSize;
		
		grids3D[g].getTemplateInfo(tmpImg, tmpPlane, tmpCenter, tmpSize);
		
		int noRows, noCols;
		noRows = grids3D[g].getNoRows();
		noCols = grids3D[g].getNoColumns();
		
		Vec3f horTrans, verTrans;
		grids3D[g].getTransformation(horTrans, verTrans);
		
		Vec3f minCoord, maxCoord;
		Vec3f xAxis, yAxis;
		Plane pln = grids3D[g].getGridPlane();
		pln.getBoundingBox(minCoord, maxCoord);
		xAxis = pln.getPlaneXAxis(); yAxis = pln.getPlaneYAxis();
		
		float xIncr = dot(xAxis, (maxCoord-minCoord)) / rectImgWidth;
		float yIncr =  dot(yAxis, (maxCoord-minCoord)) / rectImgHeight;
		
		int w = tmpSize[0]; int h = tmpSize[1];
		vector<Vec3f> templateCoordinates;
		Vec3f raydir;
		
		Vec3f p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0], tmpCenter[1], 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		
		int colIndex = -1;
		int rowIndex = -1;
		float minDist;
		for(int r=0; r<noRows; r++)
		{
			for(int c=0; c<noCols; c++)
			{
				Vec3f center = grids3D[g].getElementCenter(r,c);
				if(colIndex == -1 && rowIndex == -1)
				{
					colIndex = c;
					rowIndex = r;
					minDist = (center-p).length();
				}
				else if((center-p).length() < minDist)
				{
					colIndex = c;
					rowIndex = r;
					minDist = (center-p).length();
				}
			}
		}
		
		p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0] - w/2, tmpCenter[1] - h/2, 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		templateCoordinates.push_back(p - horTrans*colIndex - verTrans*rowIndex);
		
		p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0] + w/2, tmpCenter[1] - h/2, 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		templateCoordinates.push_back(p - horTrans*colIndex - verTrans*rowIndex);
		
		p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0] + w/2, tmpCenter[1] + h/2, 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		templateCoordinates.push_back(p - horTrans*colIndex - verTrans*rowIndex);
		
		p = rectifyingHomographies[tmpImg][tmpPlane] * Vec3f(tmpCenter[0] - w/2, tmpCenter[1] + h/2, 1.0);
		p[0] /= p[2]; p[1] /= p[2]; p[2] = 1.0;
		raydir = cameras[tmpImg].unprojectPixelwithDepth(p, 0) - shrink4To3(cameras[tmpImg].getCenter());
		p = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[tmpImg].getCenter()), raydir);
		templateCoordinates.push_back(p -horTrans*colIndex - verTrans*rowIndex);
		
		vector<vector<Line> > contourLines;
		contourLines.resize(4);
		
		PointCloud *tmpPc = new PointCloud();
		
		for(int i=0; i<noImages; i++)
		{
			if(!grids3D[g].isVisibleInImage(i))
				continue;
			
			//if(i != tmpImg)
			//	continue;
			
			//printf("******image %d*****\n", i);
			for(int r=0; r<noRows; r++)
			{
				for(int c=0; c<noCols; c++)
				{
					Vec3f v1 = templateCoordinates[0] + horTrans*c + verTrans*r;
					Vec3f v2 = templateCoordinates[1] + horTrans*c + verTrans*r;
					Vec3f v3 = templateCoordinates[2] + horTrans*c + verTrans*r;
					Vec3f v4 = templateCoordinates[3] + horTrans*c + verTrans*r;
				
					Vec3f v1Proj = cameras[i].project(expand3To4(v1), 0);
					Vec3f v2Proj = cameras[i].project(expand3To4(v2), 0);
					Vec3f v3Proj = cameras[i].project(expand3To4(v3), 0);
					Vec3f v4Proj = cameras[i].project(expand3To4(v4), 0);
					
					//side 1
					if(v1Proj[0] >= 0.0 && v1Proj[0] < imageSizes[i][0] && v1Proj[1] >= 0.0 && v1Proj[1] < imageSizes[i][1])
					{
						Line l;
						index = lineDetector->findBestCompatibleEdge(lines[i], shrink3To2(v1Proj), shrink3To2(v2Proj), score, boundaryThresh);
						if(index != -1)
						{
							l = lines[i][index];
							//printf("plot([%f %f], [%f %f], 'r');\n", l.getStartPoint2D()[0], l.getEndPoint2D()[0], l.getStartPoint2D()[1], l.getEndPoint2D()[1]);
							Vec3f dir = (cameras[i].unprojectPixelwithDepth(expand2To3(l.getStartPoint2D()), 0)) - shrink4To3(cameras[i].getCenter());
							Vec3f start3D = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[i].getCenter()), dir);
							dir = (cameras[i].unprojectPixelwithDepth(expand2To3(l.getEndPoint2D()), 0)) - shrink4To3(cameras[i].getCenter());
							Vec3f end3D = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[i].getCenter()), dir);
							tmpPc->addPoint(start3D);
							tmpPc->addPoint(end3D);
							
							start3D = start3D - horTrans*c - verTrans*r;
							end3D = end3D - horTrans*c - verTrans*r;
							
							float xStart = dot(xAxis, start3D - minCoord) / xIncr;
							float yStart = dot(yAxis, start3D - minCoord) / yIncr;
							
							float xEnd = dot(xAxis, end3D - minCoord) / xIncr;
							float yEnd = dot(yAxis, end3D - minCoord) / yIncr;
							
							contourLines[0].push_back(Line(Vec2f(xStart, yStart), Vec2f(xEnd, yEnd)));
							
							printf("plot([%f %f], [%f %f], 'r');\n", xStart, xEnd, yStart, yEnd);
						}
					}
					
					//side2
					if(v2Proj[0] >= 0.0 && v2Proj[0] < imageSizes[i][0] && v3Proj[1] >= 0.0 && v3Proj[1] < imageSizes[i][1])
					{
						Line l;
						index = lineDetector->findBestCompatibleEdge(lines[i], shrink3To2(v2Proj), shrink3To2(v3Proj), score, boundaryThresh);
						if(index != -1)
						{
							l = lines[i][index];
							//printf("plot([%f %f], [%f %f], 'r');\n", l.getStartPoint2D()[0], l.getEndPoint2D()[0], l.getStartPoint2D()[1], l.getEndPoint2D()[1]);
							Vec3f dir = (cameras[i].unprojectPixelwithDepth(expand2To3(l.getStartPoint2D()), 0)) - shrink4To3(cameras[i].getCenter());
							Vec3f start3D = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[i].getCenter()), dir);
							dir = (cameras[i].unprojectPixelwithDepth(expand2To3(l.getEndPoint2D()), 0)) - shrink4To3(cameras[i].getCenter());
							Vec3f end3D = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[i].getCenter()), dir);
							tmpPc->addPoint(start3D);
							tmpPc->addPoint(end3D);
							
							start3D = start3D - horTrans*c - verTrans*r;
							end3D = end3D - horTrans*c - verTrans*r;
							
							float xStart = dot(xAxis, start3D - minCoord) / xIncr;
							float yStart = dot(yAxis, start3D - minCoord) / yIncr;
							
							float xEnd = dot(xAxis, end3D - minCoord) / xIncr;
							float yEnd = dot(yAxis, end3D - minCoord) / yIncr;
							
							contourLines[1].push_back(Line(Vec2f(xStart, yStart), Vec2f(xEnd, yEnd)));
							
							printf("plot([%f %f], [%f %f], 'r');\n", xStart, xEnd, yStart, yEnd);
						}
					}
					
					//side3
					if(v3Proj[0] >= 0.0 && v4Proj[0] < imageSizes[i][0] && v3Proj[1] >= 0.0 && v4Proj[1] < imageSizes[i][1])
					{
						Line l;
						index = lineDetector->findBestCompatibleEdge(lines[i], shrink3To2(v3Proj), shrink3To2(v4Proj), score, boundaryThresh);
						if(index != -1)
						{
							l = lines[i][index];
							//printf("plot([%f %f], [%f %f], 'r');\n", l.getStartPoint2D()[0], l.getEndPoint2D()[0], l.getStartPoint2D()[1], l.getEndPoint2D()[1]);
							Vec3f dir = (cameras[i].unprojectPixelwithDepth(expand2To3(l.getStartPoint2D()), 0)) - shrink4To3(cameras[i].getCenter());
							Vec3f start3D = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[i].getCenter()), dir);
							dir = (cameras[i].unprojectPixelwithDepth(expand2To3(l.getEndPoint2D()), 0)) - shrink4To3(cameras[i].getCenter());
							Vec3f end3D = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[i].getCenter()), dir);
							tmpPc->addPoint(start3D);
							tmpPc->addPoint(end3D);
							
							start3D = start3D - horTrans*c - verTrans*r;
							end3D = end3D - horTrans*c - verTrans*r;
							
							float xStart = dot(xAxis, start3D - minCoord) / xIncr;
							float yStart = dot(yAxis, start3D - minCoord) / yIncr;
							
							float xEnd = dot(xAxis, end3D - minCoord) / xIncr;
							float yEnd = dot(yAxis, end3D - minCoord) / yIncr;
							
							contourLines[2].push_back(Line(Vec2f(xStart, yStart), Vec2f(xEnd, yEnd)));
							
							printf("plot([%f %f], [%f %f], 'r');\n", xStart, xEnd, yStart, yEnd);
						}
					}
					
					//side4
					if(v4Proj[0] >= 0.0 && v4Proj[0] < imageSizes[i][0] && v1Proj[1] >= 0.0 && v1Proj[1] < imageSizes[i][1])
					{
						Line l;
						index = lineDetector->findBestCompatibleEdge(lines[i], shrink3To2(v4Proj), shrink3To2(v1Proj), score, boundaryThresh);
						if(index != -1)
						{
							l = lines[i][index];
							//printf("plot([%f %f], [%f %f], 'r');\n", l.getStartPoint2D()[0], l.getEndPoint2D()[0], l.getStartPoint2D()[1], l.getEndPoint2D()[1]);
							Vec3f dir = (cameras[i].unprojectPixelwithDepth(expand2To3(l.getStartPoint2D()), 0)) - shrink4To3(cameras[i].getCenter());
							Vec3f start3D = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[i].getCenter()), dir);
							dir = (cameras[i].unprojectPixelwithDepth(expand2To3(l.getEndPoint2D()), 0)) - shrink4To3(cameras[i].getCenter());
							Vec3f end3D = grids3D[g].intersectRayWithGrid(shrink4To3(cameras[i].getCenter()), dir);
							tmpPc->addPoint(start3D);
							tmpPc->addPoint(end3D);
							
							start3D = start3D - horTrans*c - verTrans*r;
							end3D = end3D - horTrans*c - verTrans*r;
							
							float xStart = dot(xAxis, start3D - minCoord) / xIncr;
							float yStart = dot(yAxis, start3D - minCoord) / yIncr;
							
							float xEnd = dot(xAxis, end3D - minCoord) / xIncr;
							float yEnd = dot(yAxis, end3D - minCoord) / yIncr;
							
							contourLines[3].push_back(Line(Vec2f(xStart, yStart), Vec2f(xEnd, yEnd)));
							
							printf("plot([%f %f], [%f %f], 'r');\n", xStart, xEnd, yStart, yEnd);
						}
					}
				}
			}
		}
		
		tmpPc->writeMesh("/Users/ceylan/Desktop/pc1.ply");
		
		vector<Vec3f> bestLines;
		for(int i=0; i<contourLines.size(); i++)
		{
			vector<float> weights;
			weights.resize(contourLines[i].size(), 1.0);
			Vec3f bestLine;
			GeometryUtils::solveForBestFitting2DLine(contourLines[i], weights, bestLine);
			bestLines.push_back(bestLine);
		}
		
		vector<Vec2f> intersectionPts;
		GeometryUtils::findIntersectionOfContourLines(bestLines, intersectionPts);
		
		
		vector<Vec3f> corners;
		for(int i=0; i<intersectionPts.size(); i++)
		{
			printf("intersection pt:%f %f\n", intersectionPts[i][0], intersectionPts[i][1]);
			Vec3f pt = minCoord + xAxis*xIncr*intersectionPts[i][0] + yAxis*yIncr*intersectionPts[i][1];
			corners.push_back(pt);
		}
		
		/*if(g==1)
		{
			TranslationalGrid3D grid(5, 1);
			grid.setGridPlane(grids3D[g].getGridPlane());
			Vec3f h, v;
			grids3D[g].getTransformation(h, v);
			grid.setTransformations(h, v);
			grids3D[g] = grid;
			noRows = 1;
			noCols = 5;
		}*/
		
		for(int r=0; r<noRows; r++)
		{
			for(int c=0; c<noCols; c++)
			{
				vector<Vec3f> vertices;
				for(int i=0; i<corners.size(); i++)
				{
					vertices.push_back(corners[i] + horTrans*c + verTrans*r);
				}
				grids3D[g].setElement(r, c, vertices);
				
			}
		}
		
		grids3D[g].createElementMesh();
	}
	
	for(int i=0; i<noImages; i++)
	{
		Img tmpImg;
		tmpImg.resize(imageSizes[i][0], imageSizes[i][1], Color(1.0, 1.0, 1.0));
		for(int g=0; g<grids3D.size(); g++)
		{
			if(!grids3D[g].isVisibleInImage(i))
				continue;
			
			int noRows = grids3D[g].getNoRows();
			int noCols = grids3D[g].getNoColumns();
			
			for(int r=0; r<noRows; r++)
			{
				for(int c=0; c<noCols; c++)
				{
					vector<Vec3f> vertices = grids3D[g].getElement(r,c);
					for(int v=0; v<vertices.size()-1; v++)
					{
						Vec3f v1 = cameras[i].project(expand3To4(vertices[v]), 0);
						Vec3f v2 = cameras[i].project(expand3To4(vertices[v+1]),0);
						if(v1[0] >= 0.0 && v1[0]<imageSizes[i][0] && v1[1]>=0 && v1[1]<imageSizes[i][1] &&
						   v2[0] >= 0.0 && v2[0]<imageSizes[i][0] && v2[1]>=0 && v2[1]<imageSizes[i][1])
						{
							Line l(shrink3To2(v1), shrink3To2(v2));
							l.fillLineMask(imageSizes[i][0], imageSizes[i][1], &tmpImg);
						}
					}
				}
			}
		}
		
		stringstream ss;
		ss << mainDirectory << "/results/finalRepetitions_" << i << ".png";
		string filename;
		ss >> filename;
		tmpImg.write(filename);
	}
}