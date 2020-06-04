/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WienerEffect.cpp
 *
 * Description: This file contains the implementation of the
 * WienerEffect class.
 *****************************************************************************/

#include "WienerEffect.h"

namespace siena
{

/**
 * Constructor.
 */
WienerEffect::WienerEffect(
	const EffectInfo * pEffectInfo) :
		CovariateDependentContinuousEffect(pEffectInfo)
{
}


/**
 * Returns how much this effect contributes to the change in the
 * continuous behavior.
 */
double WienerEffect::calculateChangeContribution(int actor)
{
	return 0; // no contribution of the WienerEffect to the 
			  // deterministic part of the Bergstrom update
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double WienerEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;

	if (!this->missingCovariateEitherEnd(ego, this->period()))
	{
		statistic = currentValues[ego] - this->covariateValue(ego);
		statistic *= statistic;
	//	statistic = currentValues[ego] * currentValues[ego];
	}

	return statistic;
}

}
