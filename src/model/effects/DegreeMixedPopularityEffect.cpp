/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DegreeMixedPopularityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DegreeMixedPopularityEffect class.
 *****************************************************************************/

//#include <string>
//#include <stdexcept>
#include "DegreeMixedPopularityEffect.h"
#include "network/Network.h"
//#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
DegreeMixedPopularityEffect::DegreeMixedPopularityEffect(
	const EffectInfo * pEffectInfo, bool direction) :
		TwoNetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldirection = direction;
}

/**
 * Calculates the mixed popularity degree, for use in both calculateChangeContribution
 * and egoStatistic
 */
	int DegreeMixedPopularityEffect::calculateMixedPopDegree(int actor) const
{
	int statistic = 0;
	IncidentTieIterator iter;

	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();

	if (this->ldirection) // "out"
	{
		iter = pFirstNetwork->outTies(actor);
	}
	else // "in"
	{
		iter = pFirstNetwork->inTies(actor);
	}

	for ( ; iter.valid(); iter.next())
	{
		statistic+= (pSecondNetwork->inDegree(iter.actor()));
	}
	return statistic;
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double DegreeMixedPopularityEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->calculateMixedPopDegree(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double DegreeMixedPopularityEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego] * this->calculateMixedPopDegree(ego);
}

}
