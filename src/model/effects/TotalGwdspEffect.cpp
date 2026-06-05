/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalGwdspEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * TotalGwdspEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "TotalGwdspEffect.h"
#include "network/Network.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include "model/tables/ConfigurationTable.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "NetworkDependentBehaviorEffect.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{


/**
 * Constructor.
 */
TotalGwdspEffect::TotalGwdspEffect(const EffectInfo * pEffectInfo, bool forward, bool nc) :
	NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->linternalEffectParameter = pEffectInfo->internalEffectParameter();
	this->lweight = -0.01 * this->linternalEffectParameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = (1 - exp(this->lweight));
	this->lforward = forward;
	if (this->linternalEffectParameter < 0)
	{
		throw runtime_error("Gwdsp must have nonnegative internal effect parameter");
	}
	this->lnc = nc;
}



/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void TotalGwdspEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkDependentBehaviorEffect::initialize(pData, pState, period, pCache);
	// initialize the vector with weights for the GWDSP statistic
	// this is done several times during one estimation run (not elegant,
	// but not computationally burdensome either)
	double pow = 1;
	int m = this->pNetwork()->m();
	this->lcumulativeWeight.resize(m); // default values 0
	for (int i = 1; i < m; i++)
	{
		pow *= this->lexpfactor;
		this->lcumulativeWeight[i] = this->lexpmweight * (1 - pow);
	}
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TotalGwdspEffect::calculateChangeContribution(int actor,
		int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();
	
	if (pNetwork->outDegree(actor) > 0) 
	{

		this->lpInitialisedTable = 0;

		if (lforward)
					this->lpInitialisedTable = this->pTwoPathTable();
				else
					this->lpInitialisedTable = this->pInStarTable();

		double sum = 0;
		for (int j = 0; j < this->n(); j++) //inefficient?
		{
			int twoc = 0;
			if (j != actor)
			{
				twoc = this->lpInitialisedTable->get(j);
				sum += this->lcumulativeWeight[twoc];
			}
		}
		contribution = difference * sum;
	}
	return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego.
 */
double TotalGwdspEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	for (int j = 0; j < this->pNetwork()->n(); j++) // was m() until and including 1.3.11.
	{
		if (j != ego)
		{	
			int pathCount = 0;
			if (this->lforward) // tables can not be used here because ego is not preprocessed?
				pathCount = pNetwork->twoPathCount(ego, j);
			else
				pathCount = pNetwork->inTwoStarCount(ego, j);
			statistic += this->lcumulativeWeight[pathCount];
		}
	}
	statistic *= lnc ? currentValues[ego] : currentValues[ego] + this->overallCenterMean();
	return statistic;
}

}
