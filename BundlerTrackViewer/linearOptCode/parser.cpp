struct gridShift
{
	int imageIndex1;
	int imageIndex2;
	int templateId;
	int planeIndex1;
	int planeIndex2;
	int gridIndex1;
	int gridIndex2;
	int shiftCol;
	int shiftRow;
	
	gridShift()
	{
		shiftCol = 0;
		shiftRow = 0;
	}
	
	bool operator==(const gridShift &rhs) const
	{
		if(rhs.imageIndex1 == imageIndex1 && rhs.imageIndex2 == imageIndex2 &&
			rhs.planeIndex1 == planeIndex1 && rhs.planeIndex2 == planeIndex2 &&
			rhs.gridIndex1 == gridIndex1 && rhs.gridIndex2 == gridIndex2 &&
			rhs.shiftCol == shiftCol && rhs.shiftRow == shiftRow)
				return true;
		else
			return false;
	}
};

struct imageTransformation
{
	int ind1;
	int ind2;
	bool gridMatch;
	vector<gridShift> gridShifts;
	float scale;
	Vec2f translation;
	Matrix3f fundMatrix;
	float score;
	int group;
	bool used;
	bool addedToGraph;
	
	imageTransformation()
	{
		ind1 = -1;
		ind2 = -1;
		gridMatch = false;
		scale = 1.0;
		translation = Vec2f(0.0, 0.0);
		fundMatrix.setIdentity();
		group = -1;
		used = false;
		addedToGraph = false;
		gridShifts.clear();
	}
	
	bool operator<(const imageTransformation &rhs) const
	{
		if(score > rhs.score)
			return true;
		else
			return false;
	}
	
	bool operator==(const imageTransformation &rhs) const
	{
		int g = 0;
		for(; g<rhs.gridShifts.size(); g++)
		{
			if(find(gridShifts.begin(), gridShifts.end(),rhs.gridShifts[g]) == gridShifts.end())
				break;
		}
		if(g == rhs.gridShifts.size())
			return true;
		else 
			return false;
	}
};

#include <fstream>
#include "GraphWrapper.h"

void BundlerManager::convertAlignmentsToMatches(string alignmentFileName)
{
	int numImages = 0; //to be filled by you
	int noTemplates = 0; //to be filled by you
	
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
		vector<double> fundMatrix;
		fundMatrix.resize(9);
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
		fin >> tr.fundMatrix[0] >> tr.fundMatrix[1] >> tr.fundMatrix[2];
		fin >> tr.fundMatrix[3] >> tr.fundMatrix[4] >> tr.fundMatrix[5];
		fin >> tr.fundMatrix[6] >> tr.fundMatrix[7] >> tr.fundMatrix[8];
		if(var > 0.2)
			continue;
		
		gw.addEdge(img1, img2, 1.0/weight);
		
		alignments[img1][img2] = tr;
	}
	
	
	gw.findMinSpanningTree(sourceVertices, targetVertices);

	//now you can filter the feature matches based on the alignments of the min span tree
}
