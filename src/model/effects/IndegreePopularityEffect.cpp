/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreePopularityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IndegreePopularityEffect class.
 *****************************************************************************/

#include "IndegreePopularityEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "data/NetworkLongitudinalData.h"
#include "model/EffectInfo.h"
#include "data/Data.h"

namespace siena
{

/**
 * Constructor.
 */
IndegreePopularityEffect::IndegreePopularityEffect(
	const EffectInfo * pEffectInfo, bool root, bool centered): 
										NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
	this->lcentered = centered;
	this->lcentering = 0.0;
	this->lvariableName = pEffectInfo->variableName();
// centering and root cannot occur simultaneously
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void IndegreePopularityEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	if (this->lcentered)
	{
		NetworkLongitudinalData * pNetworkData =
				pData->pNetworkData(this->lvariableName);
		this->lcentering = pNetworkData->averageInDegree();
	}
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double IndegreePopularityEffect::calculateContribution(int alter) const
{
	int degree = this->pNetwork()->inDegree(alter);

	if (!this->outTieExists(alter))
	{
		// The indegree will increase after introducing the tie
		degree++;
	}

	double change = degree - this->lcentering;

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(degree);
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double IndegreePopularityEffect::tieStatistic(int alter)
{
	const Network * pNetwork = this->pNetwork();
	int degree = pNetwork->inDegree(alter);
	double statistic;

	if (this->lroot)
	{
		statistic = this->lsqrtTable->sqrt(degree);
	}
	else
	{
		statistic = degree - this->lcentering;
	}

	return statistic;
}

}
