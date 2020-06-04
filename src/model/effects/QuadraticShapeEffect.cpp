/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: QuadraticShapeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * QuadraticShapeEffect class.
 *****************************************************************************/

#include "QuadraticShapeEffect.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
QuadraticShapeEffect::QuadraticShapeEffect(const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double QuadraticShapeEffect::calculateChangeContribution(int actor,
		int difference)
{
	return (2 * this->centeredValue(actor) + difference) * difference;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double QuadraticShapeEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego] * currentValues[ego];
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double QuadraticShapeEffect::endowmentStatistic(const int * difference,
	double * currentValues)
{
	double statistic = 0;
	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
		{
			statistic += currentValues[i] * currentValues[i]
				- (currentValues[i] + difference[i]) *
				(currentValues[i] + difference[i]);
		}
	}

	return statistic;
}

}
