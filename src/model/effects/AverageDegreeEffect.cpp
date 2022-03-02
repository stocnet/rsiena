/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageDegreeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageDegreeEffect class.
 *****************************************************************************/

#include "AverageDegreeEffect.h"
#include "network/OneModeNetwork.h"
#include "data/NetworkLongitudinalData.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"
#include "data/Data.h"

namespace siena
{

/**
 * Constructor.
 */
AverageDegreeEffect::AverageDegreeEffect(
	const EffectInfo * pEffectInfo): NetworkEffect(pEffectInfo)
{
	this->lcentering = pEffectInfo->internalEffectParameter();
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void AverageDegreeEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
}



/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double AverageDegreeEffect::calculateContribution(int alter) const
{
	double change = 0;
	for (int j = 0; j < this->pNetwork()->n(); j++)
	{
		change = change + this->pNetwork()->outDegree(j);
	}
	change = (change + this->pNetwork()->outDegree(this->ego()))/this->pNetwork()->n()
				 - this->lcentering;
	if (this->outTieExists(alter))
	{
		change = change - (2/this->pNetwork()->n());
	}
	else
	{
		change = change + (2/this->pNetwork()->n());
	}
	return change;
}


/**
 * The contribution of ego to the statistic.
 */
double AverageDegreeEffect::egoStatistic(int ego, const Network * pNetwork)
{
	double totDegree = 0;
	for (int j = 0; j < pNetwork->n(); j++)
	{
		totDegree = totDegree + pNetwork->outDegree(j);
	}
	return pNetwork->outDegree(this->ego())*
			((totDegree/pNetwork->n()) - this->lcentering);
}
}
