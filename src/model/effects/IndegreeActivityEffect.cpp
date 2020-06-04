/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreeActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IndegreeActivityEffect class.
 *****************************************************************************/

#include <cmath>
#include "data/NetworkLongitudinalData.h"
#include "IndegreeActivityEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"
#include "data/Data.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
IndegreeActivityEffect::IndegreeActivityEffect(
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
void IndegreeActivityEffect::initialize(const Data * pData,
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
double IndegreeActivityEffect::calculateContribution(int alter) const
{
	double change =
		this->pNetwork()->inDegree(this->ego()) - this->lcentering;

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(change);
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double IndegreeActivityEffect::tieStatistic(int alter)
{
	const Network * pNetwork = this->pNetwork();
	int inDegree = pNetwork->inDegree(this->ego());
	double statistic = inDegree - this->lcentering;;

	if (this->lroot)
	{
		statistic = this->lsqrtTable->sqrt(inDegree);
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double IndegreeActivityEffect::endowmentStatistic(Network * pLostTieNetwork)
{
	double statistic = 0;
	const Network* pStart = this->pData()->pNetwork(this->period());
	int n = pStart->n();
	for (int i = 0; i < n; i++)
	{
		int indeg = pStart->inDegree(i);
		double rindeg = indeg - this->lcentering;
		if (this->lroot)
		{
			rindeg = this->lsqrtTable->sqrt(indeg);
		}
		statistic += rindeg * pLostTieNetwork->outDegree(i);
	}
	return statistic;
}

}
