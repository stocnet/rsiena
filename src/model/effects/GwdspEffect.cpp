/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwdspEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwdspEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "GwdspEffect.h"
#include "network/Network.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"


using namespace std;

namespace siena
{


/**
 * Constructor.
 */
GwdspEffect::GwdspEffect(const EffectInfo * pEffectInfo,
			EgocentricConfigurationTable * (NetworkCache::*pTable)() const,
			bool forward):
	NetworkEffect(pEffectInfo)
{
	this->lparameter = pEffectInfo->internalEffectParameter();
	this->lweight = -0.01 * this->lparameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = (1 - exp(this->lweight));
	this->lpTable = pTable;
	this->lforward = forward;
	if (this->lparameter < 0)
	{
		throw runtime_error("Gwdsp must have nonnegative internal effect parameter");
	}
//	this->lpNetwork = 0;
	this->lpNetworkCache = 0;
}



/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void GwdspEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	this->lpNetworkCache = pCache->pNetworkCache(this->pNetwork());
	this->lpInitialisedTable = (*this->pNetworkCache().*lpTable)();

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
double GwdspEffect::calculateContribution(int alter) const
{
	IncidentTieIterator iter;
	const Network * pNetwork = this->pNetwork();

	if (this->lforward) // "FF"
	{
		iter = pNetwork->outTies(alter);
	}
	else // "FB"
	{
		iter = pNetwork->inTies(alter);
	}

	double contribution = 0;

	for (  ; iter.valid(); iter.next())
	{
		if (iter.actor() != this->ego())
		{
			int twoc = this->lpInitialisedTable->get(iter.actor());
			if (this->outTieExists(alter))
			{
				contribution += this->lcumulativeWeight[twoc] -
						this->lcumulativeWeight[twoc-1];
			}
			else
			{
				contribution += this->lcumulativeWeight[twoc+1] -
						this->lcumulativeWeight[twoc];
			}
		}
	}
	return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego.
 */
double GwdspEffect::egoStatistic(int ego,
	const Network * pNetwork)
{
	double statistic = 0;

	for (int h = 0; h < this->pNetwork()->m(); h++)
	{
		if (h != ego)
		{
			statistic += 
					this->lcumulativeWeight[this->lpInitialisedTable->get(h)];
		}
// ego is implicit in this->lpInitialisedTable
	}

	return statistic;
}

}
