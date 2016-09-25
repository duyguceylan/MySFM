/*
 *  TextureManager.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 4/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "TextureManager.h"
#include "../MRFUtils/MRFController.h"
#include "Image.h"
#include "MathUtils.h"

vector<vector<Color> > TextureManager::colors = vector<vector<Color> > ();
vector<vector<float> > TextureManager::gradients = vector<vector<float> > ();
vector<vector<bool> > TextureManager::colorMasks = vector<vector<bool> > ();
vector<vector<int> > TextureManager::neighborRelations = vector<vector<int> > ();
vector<int> TextureManager::labelIndices = vector<int>();

float TextureManager::computePlaneSmoothnessForMRF(int pix1, int pix2, int label1, int label2)
{
	float maxCost = 5.0;
	
	if(find(neighborRelations[pix1].begin(), neighborRelations[pix1].end(), pix2) == neighborRelations[pix1].end())
	{
		printf("Error in neighboring!\n");
		return 0.0;
	}
	
	if(TextureManager::colorMasks[TextureManager::labelIndices[pix1]][label1]==false || TextureManager::colorMasks[TextureManager::labelIndices[pix2]][label2]==false ||
	   TextureManager::colorMasks[TextureManager::labelIndices[pix1]][label2]==false || TextureManager::colorMasks[TextureManager::labelIndices[pix2]][label1]==false)
	{
		//printf("s3\n");
		return maxCost;
	}
	
	if(label1 == label2)
		return 0.0;
	
	else
	{
		float diff1 = (TextureManager::colors[TextureManager::labelIndices[pix1]][label1] - TextureManager::colors[TextureManager::labelIndices[pix1]][label2]).length();
		float diff2 = (TextureManager::colors[TextureManager::labelIndices[pix2]][label1] - TextureManager::colors[TextureManager::labelIndices[pix2]][label2]).length();
		//float diff3 = abs(TextureManager::gradients[label1][pix1] - TextureManager::gradients[label2][pix1]);
		//float diff4 = abs(TextureManager::gradients[label1][pix2] - TextureManager::gradients[label2][pix2]);
		
		return diff1+diff2;
		//return diff1+diff2+diff3+diff4;
	}
	
	Color c1 = TextureManager::colors[pix1][label1];
	Color c2 = TextureManager::colors[pix2][label2];
	//printf("s2\n");
	return (c1-c2).length();
}

vector<int> TextureManager::detectAndFillOcclusions(int noLabels, int noPixels, int refIndex, int width, int height, vector<int> &pixels, vector<vector<Color> > &colors_, vector<vector<bool> > &colorMasks_, 
											 vector<int> &initialLabeling, vector<Color> &boundaryPixels, bool smoothBoundary)
{
	vector<float> dataEnergies;
	
	dataEnergies.resize(noPixels*noLabels, 0.0);
	
	TextureManager::colors.clear();
	TextureManager::colorMasks.clear();
	TextureManager::gradients.clear();
	TextureManager::colors.resize(colors_.size());
	//TextureManager::gradients.resize(noLabels);
	TextureManager::colorMasks.resize(colors_.size());
	TextureManager::neighborRelations.clear();
	TextureManager::labelIndices.clear();
	TextureManager::labelIndices.resize(noPixels);
	
	for(int i=0; i<colors_.size(); i++)
	{
		for(int j=0; j<colors_[i].size(); j++)
		{
			TextureManager::colors[i].push_back(colors_[i][j]);
			TextureManager::colorMasks[i].push_back(colorMasks_[i][j]);
		}
	}
	
	/*for(int i=0; i<noLabels; i++)
	{
		vector<float> derX;
		vector<float> derY;
		derX.resize(width*height, 0.0);
		derY.resize(width*height, 0.0);
		
		//x-derivative
		for(int y=0;y<height;y++)
		{
			for(int x=0;x<(width-1);x++)
			{
				derX[x*height+y] = (colors[x*height+y][i] - colors[(x+1)*height+y][i]).length();
			}
		}
		
		//y-derivative
		for(int x=0;x<width;x++)
		{
			int pos = x*height;
			for(int y=0; y<(height-1); y++)
			{
				derY[pos+y] = (colors[pos+y][i]-colors[pos+y+1][i]).length();
			}
		}
		
		for(int x=0; x<width; x++)
		{
			for(int y=0; y<height; y++)
			{
				if(pixels[x*height+y] != -1)
					TextureManager::gradients[i].push_back(derX[x*height+y] + derY[x*height+y]);
			}
		}
	}*/
	
	
	float maxCost = exp(3.0);
	
	//fill data cost
	#pragma omp parallel for
	for(int x=0; x<width; x++)
	{
		for(int y=0; y<height; y++)
		{
			int index = x*height+y;
			int varIndex = pixels[index];
			
			if(varIndex == -1)
				continue;
			
			TextureManager::labelIndices[varIndex] = (index);
			
			for(int i=0; i<noLabels; i++)
			{
				
				if(colorMasks_[index][i] == false)
				{
					//assign a high data cost
					dataEnergies[varIndex*noLabels+i] = maxCost;
					continue;
				}
				Color refColor = colors_[index][i];
				int count = 0;
				for(int j=0; j<noLabels; j++)
				{
					if(colorMasks_[index][j])
					{
						dataEnergies[varIndex*noLabels+i] += (refColor - colors_[index][j]).length();
						count += 1;
					}
				}
				dataEnergies[varIndex*noLabels+i] = dataEnergies[varIndex*noLabels+i] / (float)count;
				
				if(smoothBoundary)
				{
					float dist  = 0.0;
					//top left corner
					if(x==0 && y==0)
					{
						//5 neighbors
						dist += (refColor-boundaryPixels[0]).length();
						dist += (refColor-boundaryPixels[1]).length();
						dist += (refColor-boundaryPixels[2]).length();
						dist += (refColor-boundaryPixels[width+2]).length();
						dist += (refColor-boundaryPixels[width+2+2]).length();
						dist /= 5.0;
						printf("dist1:%f\n", dist);
					}
					//top right corner
					else if(x==width-1 && y==0)
					{
						//5 neighbors
						dist += (refColor-boundaryPixels[width-1]).length();
						dist += (refColor-boundaryPixels[width]).length();
						dist += (refColor-boundaryPixels[width+1]).length();
						dist += (refColor-boundaryPixels[width+2+1]).length();
						dist += (refColor-boundaryPixels[width+2+3]).length();
						dist /= 5.0;
						printf("dist2:%f\n", dist);
					}
					//bottom left corner
					else if(x==0 && y==height-1)
					{
						//5 neighbors
						dist += (refColor-boundaryPixels[width+2+(height-2)*2]).length();
						dist += (refColor-boundaryPixels[width+2+(height-1)*2]).length();
						dist += (refColor-boundaryPixels[width+2+height*2]).length();
						dist += (refColor-boundaryPixels[width+2+height*2+1]).length();
						dist += (refColor-boundaryPixels[width+2+height*2+2]).length();
						dist /= 5.0;
						printf("dist3:%f\n", dist);
					}
					//bottom right corner
					else if(x==width-1 && y==height-1)
					{
						dist += (refColor-boundaryPixels[width+2+(height-2)*2+1]).length();
						dist += (refColor-boundaryPixels[width+2+(height-1)*2+1]).length();
						dist += (refColor-boundaryPixels[width+2+height*2+width+1]).length();
						dist += (refColor-boundaryPixels[width+2+height*2+width]).length();
						dist += (refColor-boundaryPixels[width+2+height*2+width-1]).length();
						dist /= 5.0;
						printf("dist4:%f\n", dist);
					}
					else if(y==0)
					{
						dist += (refColor-boundaryPixels[x]).length();
						dist += (refColor-boundaryPixels[x+1]).length();
						dist += (refColor-boundaryPixels[x+2]).length();
						dist /= 3.0;
						printf("dist5:%f\n", dist);
					}
					else if(y==height-1)
					{
						dist += (refColor-boundaryPixels[width+2+height*2 + x]).length();
						dist += (refColor-boundaryPixels[width+2+height*2 + x+1]).length();
						dist += (refColor-boundaryPixels[width+2+height*2 + x+2]).length();
						dist /= 3.0;
						printf("dist6:%f\n", dist);
					}
					else if(x==0)
					{
						dist += (refColor-boundaryPixels[width+2+(y-1)*2]).length();
						dist += (refColor-boundaryPixels[width+2+y*2]).length();
						dist += (refColor-boundaryPixels[width+2+(y+1)*2]).length();
						dist /= 3.0;
						printf("dist7:%f\n", dist);
					}
					else if(x==width-1)
					{
						dist += (refColor-boundaryPixels[width+2+(y-1)*2+1]).length();
						dist += (refColor-boundaryPixels[width+2+y*2+1]).length();
						dist += (refColor-boundaryPixels[width+2+(y+1)*2+1]).length();
						dist /= 3.0;
						printf("dist8:%f\n", dist);
					}
					
					dataEnergies[varIndex*noLabels+i] = dataEnergies[varIndex*noLabels+i]+dist*3;
				}
				
				dataEnergies[varIndex*noLabels+i] = exp(dataEnergies[varIndex*noLabels+i]);
				/*if(i != refIndex)
				{
					dataEnergies[varIndex*noImages+i] += 0.1;
				}*/
			}
		}
	}
	
	//fill neighboring relations
	for(int x=0; x<width; x++)
	{
		for(int y=0; y<height; y++)
		{
			if(pixels[x*height+y] != -1)
			{
				vector<int> neighbors;
				
				if(x-1>0)
				{
					if(pixels[(x-1)*height+y] != -1)
					{
						neighbors.push_back(pixels[(x-1)*height+y]);
					}
				}
				if(x+1<width)
				{
					if(pixels[(x+1)*height+y] != -1)
					{
						neighbors.push_back(pixels[(x+1)*height+y]);
					}
				}
				if(y-1>0)
				{
					if(pixels[x*height+(y-1)] != -1)
					{
						neighbors.push_back(pixels[x*height+(y-1)]);
					}
				}
				if(y+1<height)
				{
					if(pixels[x*height+(y+1)] != -1)
					{
						neighbors.push_back(pixels[x*height+(y+1)]);
					}
				}
				
				TextureManager::neighborRelations.push_back(neighbors);
			}
		}
	}
	
	MRFController *mrfController = new MRFController();
	MRFParameters param;
    param.lambda     = 0.25;
    param.smoothExp = 2;
    param.smoothMax = 5.0;
    param.iterMax   = 10;
    param.verbose    = 1;
	param.optType	 = GraphcutExpansion;
	mrfController->setParam(param);
	
	vector<int> finalLabels;
	mrfController->performMRFWithGraphCutExpansion(dataEnergies, TextureManager::neighborRelations, TextureManager::computePlaneSmoothnessForMRF, initialLabeling, finalLabels, noPixels, noLabels);
    //mrfController->performMRFWithGraphCutExpansion(dataEnergies, vCues, hCues,finalLabels, width, height, noImages);
	delete mrfController;
	
	Img assignments;
	assignments.resize(width, height, Color(0.0, 0.0, 0.0));
	
	for(int x=0; x<width; x++)
	{
		for(int y=0; y<height; y++)
		{
			int index = x*height+y;
			int varIndex = pixels[index];
			if(varIndex == -1)
				continue;
			
			if(finalLabels[varIndex] < 0 || finalLabels[varIndex] > noLabels)
			{
				printf("Error:Invalid label assignment!\n");
			}
			else if(finalLabels[varIndex] == refIndex)
			{
				assignments.setColor(x, y, Color(0.0, 0.0, 0.0));
			}
			else
			{
				Vec3uc c = MathUtils::generateColorFromValue(finalLabels[varIndex], 0.0, noLabels);
				assignments.setColor(x, y, Color(c[0]/255.0, c[1]/255.0, c[2]/255.0));
			}
		}
	}
	assignments.write("/Users/ceylan/Desktop/assignments.png");
	
	TextureManager::colors.clear();
	TextureManager::gradients.clear();
	TextureManager::colorMasks.clear();
	TextureManager::neighborRelations.clear();
	
	return finalLabels;
}