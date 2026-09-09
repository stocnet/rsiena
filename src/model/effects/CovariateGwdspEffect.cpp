/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateGwdspEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateGwdspEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "CovariateGwdspEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
CovariateGwdspEffect::CovariateGwdspEffect(const EffectInfo * pEffectInfo,
		EgocentricConfigurationTable * (NetworkCache::*pTable)() const,
		bool forward, bool nc) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lparameter = pEffectInfo->internalEffectParameter();
	this->lweight = -0.01 * this->lparameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = (1 - exp(this->lweight));
	this->lpTable = pTable;
	this->lforward = forward;
	this->lnc = nc;
	if (this->lparameter < 0)
	{
		throw runtime_error(
			"Gwdsp must have nonnegative internal effect parameter");
	}
	this->lpNetworkCache = 0;
}


/**
 * Initializes this effect.
 */
void CovariateGwdspEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);
	this->lpNetworkCache = pCache->pNetworkCache(this->pNetwork());
	this->lpInitialisedTable = (*this->lpNetworkCache.*lpTable)();

	// Precompute the GW weights, indexed by the shared-partner count S_ih.
	// The covariate weight z_h multiplies the transform but does not change
	// the index, so S_ih <= m() as in the structural GwdspEffect.
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
 * For potential tie ego -> alter this is
 *   sum_{h in in/outTies(alter), h != ego} z_h * ( w(S_{ego,h}+/-1) - w(S_{ego,h}) ),
 * i.e. the covariate-weighted version of GwdspEffect::calculateContribution.
 */
double CovariateGwdspEffect::calculateContribution(int alter) const
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

	for ( ; iter.valid(); iter.next())
	{
		if (iter.actor() != this->ego())
		{
			int twoc = this->lpInitialisedTable->get(iter.actor());
			double increment;
			if (this->outTieExists(alter))
			{
				increment = this->lcumulativeWeight[twoc] -
						this->lcumulativeWeight[twoc - 1];
			}
			else
			{
				increment = this->lcumulativeWeight[twoc + 1] -
						this->lcumulativeWeight[twoc];
			}
			contribution += this->resolvedValue(iter.actor(), this->lnc) * increment;
		}
	}
	return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego:
 *   sum_{h != ego} z_h * w(S_{ego,h}).
 */
double CovariateGwdspEffect::egoStatistic(int ego,
	const Network * pSummationTieNetwork)
{
	double statistic = 0;

	for (int h = 0; h < this->pNetwork()->n(); h++)
	{
		if (h != ego && !this->missing(h))
		{
			double covVal = this->resolvedValue(h, this->lnc);
			statistic += covVal *
						this->lcumulativeWeight[this->lpInitialisedTable->get(h)];
		}
		// ego is implicit in this->lpInitialisedTable
	}

	return statistic;
}

}
