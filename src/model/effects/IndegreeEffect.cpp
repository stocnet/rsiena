/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IndegreeEffect class.
 *****************************************************************************/

#include <cmath>
#include "IndegreeEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
IndegreeEffect::IndegreeEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double IndegreeEffect::calculateChangeContribution(int actor,
	int difference)
{
	// The formula for the effect:
	// s_i(x) = v_i * indegree of i.
	// We need to calculate the change delta in s_i(x), if we changed
	// v_i to v_i + d (d being the given amount of change in v_i).
	// This is  d * indegree of i. This is what is calculated below.

	return difference * this->pNetwork()->inDegree(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double IndegreeEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego] * this->pNetwork()->inDegree(ego);
}


/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double IndegreeEffect::egoEndowmentStatistic(int i, const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if (difference[i] > 0)
	{
		//	statistic += currentValues[i] * this->pNetwork()->inDegree(i);
		statistic -= difference[i] * this->pNetwork()->inDegree(i);
	}

	return statistic;
}

}
