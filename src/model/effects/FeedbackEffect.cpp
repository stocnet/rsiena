/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: FeedbackEffect.cpp
 *
 * Description: This file contains the implementation of the
 * FeedbackEffect class.
 *****************************************************************************/

#include "FeedbackEffect.h"

namespace siena
{

/**
 * Constructor.
 */
FeedbackEffect::FeedbackEffect(
	const EffectInfo * pEffectInfo) :
		CovariateDependentContinuousEffect(pEffectInfo)
{
}


/**
 * Returns how much this effect contributes to the change in the
 * continuous behavior.
 */
double FeedbackEffect::calculateChangeContribution(int actor)
{
	return this->covariateValue(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double FeedbackEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;

	if (!this->missingCovariateEitherEnd(ego, this->period()))
	{
		statistic = currentValues[ego] * this->covariateValue(ego);
	}

	return statistic;
}

}
