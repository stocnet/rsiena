/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalGwdspAlterEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwdspEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "TotalGwdspAlterEffect.h"
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

// would probably be nice if GwespFunction could be used -> implement as network method?

using namespace std;

namespace siena
{


/**
 * Constructor.
 */
TotalGwdspAlterEffect::TotalGwdspAlterEffect(const EffectInfo * pEffectInfo, bool forward) :
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
//	this->lpNetwork = 0;
//	this->lpNetworkCache = 0;
}



/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void TotalGwdspAlterEffect::initialize(const Data * pData,
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
double TotalGwdspAlterEffect::calculateChangeContribution(int actor,
		int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();
	
	if (pNetwork->outDegree(actor) > 0) 
	{
		// The formula for the effect:
		// tbd

		this->lpInitialisedTable = 0;

		if (lforward)
					this->lpInitialisedTable = this->pTwoPathTable();
				else
					this->lpInitialisedTable = this->pInStarTable();

		double sumAlterValue = 0;
		// double denom = 0;
		for (int j = 0; j < this->n(); j++) //inefficient?
		{
			double alterValue = 0;
			int twoc = 0;
			if (j != actor)
			{
				if (this->lforward)
					twoc = this->lpInitialisedTable->get(j);
				else 
					twoc = this->lpInitialisedTable->get(j);
				alterValue = this->centeredValue(j) * this->lcumulativeWeight[twoc];
				// int tieValue =  this->pNetwork()->tieValue(actor, j);
				// if (((pNetwork->inDegree(j) - tieValue)> 0) && (this->ldivide2))
				// {
				// 	alterValue /= (pNetwork->inDegree(j) - tieValue);
				// }
				sumAlterValue += alterValue;
			}
		}
		contribution = difference * sumAlterValue;
		// if (denom != 0)
		// {
		// 	contribution /= denom; //what happens if denom == 0 but sumaAlterValue != 0 ?
		// }
	// 	if (this->ldivide1)
	// 	{
	// 		contribution /= pNetwork->outDegree(actor);
	// 	}
	}
	return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego.
 */
double TotalGwdspAlterEffect::egoStatistic(int ego, double * currentValues)
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
			statistic += (currentValues[j]) *
					this->lcumulativeWeight[pathCount];
		}
	}
	statistic *= (currentValues[ego]);

	return statistic;
}

}
