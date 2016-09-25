/*
 *  MRFController.cpp
 *  SymmMVS
 *
 *  Created by Duygu Ceylan on 5/31/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "MRFController.h"
#include <cmath>

#define BIGNUM 3.3e33

double round(double number)
{
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

MRFController::MRFController()
{
    initialize();
}

MRFController::MRFController(const MRFParameters& param)
{
    initialize();
    setParam(param);
}

MRFController::~MRFController()
{
    release();
}

void MRFController::initialize()
{
    // set default parameters.
    mrfParams.lambda       = 0.1f;
    mrfParams.smoothExp   = 2.0f;
    mrfParams.smoothMax   = BIGNUM;
    mrfParams.iterMax     = 5;
    mrfParams.verbose      = 1;
}

void MRFController::release() {}

void MRFController::setParam(const MRFParameters& mrfParams_)
{
    mrfParams = mrfParams_;
}

void MRFController::performWTA(vector<float>& totalEnergies, vector<float> &initialLabels, vector<float>& finalLabels, int w, int h, int noLabels)
{
	finalLabels.resize(w*h);
	for(int x=0; x<w; x++)
	{
		for(int y=0; y<h; y++)
		{
			float minEnergy = totalEnergies[(y*w+x)*noLabels];
			float minLabel = 0;
			for(int i=1; i<noLabels; i++)
			{
				if(totalEnergies[(y*w+x)*noLabels+i] < minEnergy)
				{
					minEnergy = totalEnergies[(y*w+x)*noLabels+i];
					minLabel = i;
				}
			}
			finalLabels[y*w+x] = minLabel;
		}
	}
}

/*void MRFController::performMRFWithGraphCutExpansion(vector<float>& totalEnergies, vector<float> &labels, int w, int h, int noLabels)
{
	//to access the cost of pixel (x,y) for label l: totalEnergies[(y*width+x)*noLabels + l]
	SmoothnessCost smoothCost(mrfParams.smoothExp, mrfParams.smoothMax, mrfParams.lambda);
    DataCost dataCost(&(totalEnergies[0]));
    EnergyFunction energy(&dataCost, &smoothCost);
    Expansion mrf(w, h, noLabels, &energy);
	
    vector<float> initialLabels;
	initialLabels.resize(w*h, 0.0);
    labels.resize(w*h);
	
    optimize(mrf, initialLabels, labels);
}

void MRFController::performMRFWithGraphCutExpansion(vector<float>& totalEnergies, vector<float> &vcue, vector<float> &hcue, 
													vector<int> &finalLabels, int w, int h, int noLabels)
{
	GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(w, h, noLabels);
	gc->setVerbosity(1);
	
	//data cost
	gc->setDataCost(&(totalEnergies[0]));
	
    SmoothnessCost smoothCost(mrfParams.smoothExp, mrfParams.smoothMax, mrfParams.lambda, &(hcue[0]), &(vcue[0]));
    DataCost dataCost(&(totalEnergies[0]));
    EnergyFunction energy(&dataCost, &smoothCost);
    Expansion mrf(w, h, noLabels, &energy);
	vector<int> initialLabels;
	initialLabels.resize(w*h, 0);
    finalLabels.resize(w*h);
	
    optimize(mrf, initialLabels, finalLabels);
}*/

void MRFController::performMRFWithGraphCutExpansion(vector<float>& totalEnergies, float (*costFn)(int s1, int s2, int l1, int l2),
									   vector<int> &initialLabels, vector<int> &finalLabels, int w, int h, int noLabels)
{
	GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(w, h, noLabels);
	gc->setVerbosity(1);
	
	//data cost
	gc->setDataCost(&(totalEnergies[0]));
	// smoothness comes from function pointer
	gc->setSmoothCost(costFn);
	
	/*for(int x=0; x<w; x++)
	{
		for(int y=0; y<h; y++)
		{
			gc->setLabel(y*w+x, initialLabels[y*w+x]);
		}
	}*/
	
	gc->expansion(10);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
	
	finalLabels.resize(w*h);
	for(int x=0; x<w; x++)
	{
		for(int y=0; y<h; y++)
		{
			finalLabels[y*w+x] = gc->whatLabel(y*w+x);
		}
	}
	
	delete gc;
}

void MRFController::performMRFWithGraphCutExpansion(vector<float>& totalEnergies, vector<vector<int> > &neighbors, float (*costFn)(int s1, int s2, int l1, int l2),
													vector<int> &initialLabels, vector<int> &finalLabels, int noSites, int noLabels)
{
	GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(noSites,noLabels);
	gc->setVerbosity(1);
	
	//data cost
	gc->setDataCost(&(totalEnergies[0]));
	// smoothness comes from function pointer
	gc->setSmoothCost(costFn);
	
	for(int x=0; x<noSites; x++)
	{
		gc->setLabel(x, initialLabels[x]);
	}
	
	for(int x=0; x<neighbors.size(); x++)
	{
		for(int y=0; y<neighbors[x].size(); y++)
		{
			if(x > neighbors[x][y])
			{
				gc->setNeighbors(x, neighbors[x][y], 1.0);
			}
		}
	}
	
	gc->expansion(10);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
	
	finalLabels.resize(noSites);
	for(int x=0; x<noSites; x++)
	{
		finalLabels[x] = gc->whatLabel(x);
		
	}
	
	delete gc;
}

/*void MRFController::performMRFWithICM(vector<float>& totalEnergies, vector<float> &vcue, vector<float> &hcue, 
									  vector<float> &initialLabels, vector<float> &finalLabels, int w, int h, int noLabels)
{
    SmoothnessCost smoothCost(mrfParams.smoothExp, mrfParams.smoothMax, mrfParams.lambda, &(hcue[0]), &(vcue[0]));
    DataCost dataCost(&(totalEnergies[0]));
    EnergyFunction energy(&dataCost, &smoothCost);
    ICM mrf(w, h, noLabels, &energy);
	
    finalLabels.resize(w*h);
	
    optimize(mrf, initialLabels, finalLabels);
}

void MRFController::performMRFWithTRWS(vector<float>& totalEnergies, vector<float> &vcue, vector<float> &hcue, 
									   vector<float> &initialLabels, vector<float> &finalLabels, int w, int h, int noLabels)
{
    SmoothnessCost smoothCost(mrfParams.smoothExp, mrfParams.smoothMax, mrfParams.lambda, &(hcue[0]), &(vcue[0]));
    DataCost dataCost(&(totalEnergies[0]));
    EnergyFunction energy(&dataCost, &smoothCost);
    TRWS mrf(w, h, noLabels, &energy);
	
    finalLabels.resize(w*h);
	
    optimize(mrf, initialLabels, finalLabels);
}

void MRFController::performMRFWithTRWS(vector<float>& totalEnergies, MRF::SmoothCostGeneralFn costFn, 
									   vector<int> &finalLabels, int w, int h, int noLabels)
{
	SmoothnessCost smoothCost(costFn);
    DataCost dataCost(&(totalEnergies[0]));
    EnergyFunction energy(&dataCost, &smoothCost);
    TRWS mrf(w, h, noLabels, &energy);
	
    finalLabels.resize(w*h);
	vector<int> initialLabels;
	initialLabels.resize(w*h, 0);
    optimize(mrf, initialLabels, finalLabels);
}

void MRFController::optimize(MRF& mrf, vector<int> &initialLabels, vector<int> &finalLabels)
{	
    mrf.initialize();
    mrf.clearAnswer();
	
    setLabel(mrf, initialLabels);
	
	if(mrfParams.optType == GraphcutExpansion || mrfParams.optType == GraphcutSwap)
	{
		bool randomLabelOrder[1] = {0};
		mrf.setParameters(1, &randomLabelOrder);
	}
	
    float Ed, Es, E, Eold;
    float totalTime = 0.0f;
    Ed = mrf.dataEnergy();
    Es = mrf.smoothnessEnergy();
    E  = Ed + Es;
    Eold = E;
	
	if (mrfParams.verbose) 
	{
		printf("Energy = %f (Ed=%f, Es=%f) at start (%.3f secs init)\n", E, Ed, Es, totalTime);
    }
	
	int minEnergyIter = -1;
	float minEnergy;
	
    for (int iter = 0; iter < mrfParams.iterMax; iter++) 
	{
		float t;
		mrf.optimize(1, t);
		totalTime += t;
		
		Ed = mrf.dataEnergy();
		Es = mrf.smoothnessEnergy();
		E  = Ed + Es;
		
		if(minEnergyIter == -1)
		{
			minEnergyIter = iter;
			minEnergy = E;
			getLabel(mrf, finalLabels);
		}
		else if(E < minEnergy)
		{
			minEnergy = E;
			minEnergyIter = iter;
			getLabel(mrf, finalLabels);
		}
		
		if (mrfParams.verbose) 
		{
			printf("Energy = %f (Ed=%f, Es=%f) at start (%.3f secs init)\n", E, Ed, Es, totalTime);
		}
		
		if (E == Eold) 
		{
			break;
		} 
		else if (E > Eold) 
		{
			if (mrfParams.verbose) 
			{
				printf("Warning: energy is increasing!\n");
			}
		}
	    
		Eold = E;
    }
}

void MRFController::getLabel(MRF& mrf, vector<int> &label)
{
    for (int i = 0; i < label.size(); i++) 
	{
		label[i] = float(mrf.getLabel(i));
    }
}

void MRFController::setLabel(MRF& mrf, vector<int> &label)
{
    for (int i = 0; i < label.size(); i++) {
		mrf.setLabel(i, label[i]);
    }
}*/