/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: LinearShapeNCEffect.cpp
 *
 * Description: This file contains the implementation of the
 * LinearShapeEffect class.
 *****************************************************************************/

#include "LinearShapeNCEffect.h"
#include "model/variables/BehaviorVariable.h"
#include "data/BehaviorLongitudinalData.h"

namespace siena
{

/**
 * Constructor.
 */
LinearShapeNCEffect::LinearShapeNCEffect(const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double LinearShapeNCEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double LinearShapeNCEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego] + this->overallCenterMean();
}


/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double LinearShapeNCEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if (difference[ego] > 0)
	{
		//	statistic = currentValues[ego] ;
		statistic -= difference[ego];
	}

	return statistic;
}

}
