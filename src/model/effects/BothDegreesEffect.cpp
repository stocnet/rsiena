/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BothDegreesEffect.cpp
 *
 * Description: This file contains the implementation of the
 * BothDegreesEffect class.
 *****************************************************************************/

#include <stdexcept>
#include "BothDegreesEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"
#include "Effect.h"
#include "data/Data.h"

namespace siena
{

/**
 * Constructor.
 */
BothDegreesEffect::BothDegreesEffect(
	const EffectInfo * pEffectInfo, bool centered) : NetworkEffect(pEffectInfo)
{
	this->lroot = (pEffectInfo->internalEffectParameter() >= 2);
	this->lsqrtTable = SqrtTable::instance();
	this->lcentered = centered;
	this->lcentering = 0.0;
	this->lvariableName = pEffectInfo->variableName();
	if ((this->lcentered) & (this->lroot))
	{
		throw std::logic_error(
"centering and square root may not be combined for degree activity plus popularity effect.");
	}
// centering and root cannot occur simultaneously
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void BothDegreesEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	if (this->lcentered)
	{
		NetworkLongitudinalData * pNetworkData =
				pData->pNetworkData(this->lvariableName);
		this->lcentering = 2*(pNetworkData->averageOutDegree());
	}
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double BothDegreesEffect::calculateContribution(int alter) const
{
	// see IndegreePopularityEffect, OutdegreeActivityEffect, OutdegreeActivitySqrtEffect
	// Current in-degree
	int degree1 = this->pNetwork()->inDegree(alter);

	double change2 = 0;
	// Current out-degree
	int d =	this->pNetwork()->outDegree(this->ego());
	if (this->lroot)
	{
		if (this->outTieExists(alter))
		{
			change2 =  (d * this->lsqrtTable->sqrt(d)) -
				((d-1) * this->lsqrtTable->sqrt(d-1));
		}
		else
		{
			change2 =  ((d+1) * this->lsqrtTable->sqrt(d+1)) -
				(d * this->lsqrtTable->sqrt(d));
			degree1++;
		}
	}
	else
	{
		if (this->outTieExists(alter))
		{
			change2 = 2 * d - 1 - this->lcentering;
		}
		else
		{
			change2 = 2 * d + 1 - this->lcentering;
			degree1++;
		}
	}

	double change1 = 0;
	if (this->lroot)
	{
		change1 = this->lsqrtTable->sqrt(degree1);
	}
	else
	{
		change1 = degree1;
	}

	return (change1 + change2);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double BothDegreesEffect::tieStatistic(int alter)
{
	const Network * pNetwork = this->pNetwork();
	int degree1 = pNetwork->inDegree(alter);
	int degree2 = pNetwork->outDegree(this->ego());
	double statistic;

	if (this->lroot)
	{
		statistic =  this->lsqrtTable->sqrt(degree1) +
						this->lsqrtTable->sqrt(degree2);
	}
	else
	{
		statistic = degree1 + degree2 - this->lcentering;
	}

	return statistic;
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double BothDegreesEffect::endowmentStatistic(Network * pLostTieNetwork)
{
	double statistic = 0;
	const Network* pStart = this->pData()->pNetwork(this->period());
	int n = pStart->n();
	for (int i = 0; i < n; i++)
	{
		int indeg = pStart->inDegree(i);
		double rindeg = 2*indeg;
		if (this->lroot)
		{
			rindeg = this->lsqrtTable->sqrt(indeg);
		}
		else
		{
			rindeg -= this->lcentering;
		}
		statistic += rindeg * pLostTieNetwork->outDegree(i);
	}
	return statistic;
}


}
