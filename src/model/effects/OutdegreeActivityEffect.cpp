/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreeActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutdegreeActivityEffect class.
 *****************************************************************************/

#include "OutdegreeActivityEffect.h"
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
OutdegreeActivityEffect::OutdegreeActivityEffect(
	const EffectInfo * pEffectInfo, bool centered): NetworkEffect(pEffectInfo)
{
	this->lcentered = centered;
	this->lcentering = 0.0;
	this->lvariableName = pEffectInfo->variableName();
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void OutdegreeActivityEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	if (this->lcentered)
	{
		NetworkLongitudinalData * pNetworkData =
				pData->pNetworkData(this->lvariableName);
		this->lcentering = pNetworkData->averageOutDegree();
	}
}



/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutdegreeActivityEffect::calculateContribution(int alter) const
{
	double change = 0;

	// Current out-degree
	int d = this->pNetwork()->outDegree(this->ego());

	if (this->outTieExists(alter))
	{
		// After a tie withdrawal, the new out-degree would be d-1, and
		// the new effect statistic s_i'=(d-1)^2. The current statistic
		// is s_i=d^2, so the change would be s_i-s_i' = 2d-1.
		change = 2*d - 1 - this->lcentering;
	}
	else
	{
		// When introducing a new tie, the new out-degree would be d+1, and
		// the new effect statistic s_i'=(d+1)^2. The current statistic
		// is s_i=d^2, so the change would be s_i'-s_i = 2d+1.
		change = 2*d + 1 - this->lcentering;
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double OutdegreeActivityEffect::tieStatistic(int alter)
{
	return (this->pNetwork()->outDegree(this->ego()) - this->lcentering);
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double OutdegreeActivityEffect::endowmentStatistic(Network * pLostTieNetwork)
{
	double statistic = 0;
	const Network* pStart = this->pData()->pNetwork(this->period());
	int n = pStart->n();
	for (int i = 0; i < n; i++)
	{
		int outdeg = (pStart->outDegree(i) - this->lcentering);
		statistic += outdeg * pLostTieNetwork->outDegree(i);
	}
	return statistic;
}

}
