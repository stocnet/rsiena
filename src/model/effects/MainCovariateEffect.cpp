/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MainCovariateEffect.cpp
 *
 * Description: This file contains the implementation of the
 * MainCovariateEffect class.
 *****************************************************************************/

#include "MainCovariateEffect.h"

namespace siena
{

/**
 * Constructor.
 */
MainCovariateEffect::MainCovariateEffect(
	const EffectInfo * pEffectInfo) :
		CovariateDependentBehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double MainCovariateEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->covariateValue(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double MainCovariateEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;

	if (!this->missingCovariateEitherEnd(ego, this->period()))
	{
		statistic = currentValues[ego] * this->covariateValue(ego);
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to the given ego as part of
*  the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double MainCovariateEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if (difference[ego] > 0)
	{
		statistic = - difference[ego] * this->covariateValue(ego);
	}

	return statistic;
}

}
