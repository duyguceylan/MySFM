/*
 *  MRFController.h
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/31/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _MRF_CONTROLLER_
#define _MRF_CONTROLLER_

//#include "mrf.h"
//#include "GCoptimization.h"
//#include "ICM.h"
//#include "TRW-S.h"
#include "gco/GCoptimization.h"
#include <vector>

using namespace std;

typedef enum MRFOptimizationAlgorithm
{
	IterCondModels = 0,
	GraphcutExpansion,
	GraphcutSwap,
	MaxProductBeliefPropogation,
	TreeReweightedMsgPassing,
	BPS
};

typedef struct MRFParameters
{
    float lambda;
	float smoothExp;
    float smoothMax;
    int   iterMax;
    int   verbose;
	MRFOptimizationAlgorithm optType;
};

class MRFController
{
public:
    MRFController();
    MRFController(const MRFParameters& param);
    ~MRFController();
	
	void initialize();
    void release();
	
    void setParam(const MRFParameters& param);
	
	void performWTA(vector<float>& totalEnergies, vector<float> &initialLabels, vector<float>& finalLabels, int w, int h, int noLabels);
	
    //void performMRFWithGraphCutExpansion(vector<float>& totalEnergies, vector<float> &labels, int w, int h, int noLabels);
    //void performMRFWithGraphCutExpansion(vector<float>& totalEnergies, vector<float> &vcue, vector<float> &hcue, 
	//									 vector<int> &finalLabels, int w, int h, int noLabels);
	void performMRFWithGraphCutExpansion(vector<float>& totalEnergies, float (*costFn)(int s1, int s2, int l1, int l2),
										 vector<int> &initialLabels, vector<int> &finalLabels, int w, int h, int noLabels);
	void performMRFWithGraphCutExpansion(vector<float>& totalEnergies, vector<vector<int> > &neighbors, float (*costFn)(int s1, int s2, int l1, int l2),
										 vector<int> &initialLabels, vector<int> &finalLabels, int noSites, int noLabels);
	//void performMRFWithICM(vector<float>& totalEnergies, vector<float> &vcue, vector<float> &hcue, 
	//					   vector<float> &initialLabels, vector<float> &finalLabels, int w, int h, int noLabels);
	
	//void performMRFWithTRWS(vector<float>& totalEnergies, vector<float> &vcue, vector<float> &hcue, 
	//						vector<float> &initialLabels, vector<float> &finalLabels, int w, int h, int noLabels);
	//void performMRFWithTRWS(vector<float>& totalEnergies, MRF::SmoothCostGeneralFn costFn, 
	//						vector<int> &finalLabels, int w, int h, int noLabels);
	
    //void optimize(MRF& mrf, vector<int> &initialLabels, vector<int> &finalLabels);
    
	//void getLabel(MRF& mrf, vector<int> &label);
    //void setLabel(MRF& mrf, vector<int> &label);
	
protected:
    MRFParameters mrfParams;
};

#endif