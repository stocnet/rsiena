/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocalDegreeBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ReciprocalDegreeBehaviorEffect class.
 *****************************************************************************/

#include <string>
#include <stdexcept>
#include "ReciprocalDegreeBehaviorEffect.h"
#include "network/OneModeNetwork.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
ReciprocalDegreeBehaviorEffect::ReciprocalDegreeBehaviorEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double ReciprocalDegreeBehaviorEffect::calculateChangeContribution(int actor,
	int difference)
{
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreeBehaviorEffect");
	}

	return difference * pNetwork->reciprocalDegree(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double ReciprocalDegreeBehaviorEffect::egoStatistic(int ego, double * currentValues)
{
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreeBehaviorEffect");
	}

	return currentValues[ego] * pNetwork->reciprocalDegree(ego);
}

}
