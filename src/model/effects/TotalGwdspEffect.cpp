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
#include "network/IncidentTieIterator.h"
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
	this->lparameter = pEffectInfo->internalEffectParameter();
	this->lweight = -0.01 * this->lparameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = (1 - exp(this->lweight));
	this->lforward = forward;
	if (this->lparameter < 0)
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

	this->lTwoPathCount.assign(this->pNetwork()->n(), 0);
	this->lTouched.clear();
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
	const Network * pNetwork = this->pNetwork();
	// Scatter: accumulate pathCount(ego, j) for every reachable focal actor j, by
	// walking ego's gateways h (out-ties) and, through each, the j reached (FF: out-
	// ties of h; FB: in-ties of h).
	for (IncidentTieIterator iterH = pNetwork->outTies(ego); iterH.valid(); iterH.next())
	{
		int h = iterH.actor();
		IncidentTieIterator iterJ = this->lforward ?
			pNetwork->outTies(h) : pNetwork->inTies(h);
		for (; iterJ.valid(); iterJ.next())
		{
			int j = iterJ.actor();
			if (j == ego) continue;
			if (this->lTwoPathCount[j] == 0) this->lTouched.push_back(j);
			this->lTwoPathCount[j]++;
		}
	}
	double statistic = 0;
	for (int j : this->lTouched)
	{
		statistic += this->lcumulativeWeight[this->lTwoPathCount[j]];
		this->lTwoPathCount[j] = 0;
	}
	this->lTouched.clear();
	statistic *= lnc ? currentValues[ego] + this->overallCenterMean() :
					   currentValues[ego];
	return statistic;
}

}
