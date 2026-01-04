/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: QuadraticShapeNCEffect.cpp
 *
 * Description: This file contains the implementation of the
 * QuadraticShapeNCEffect class.
 *****************************************************************************/

#include "QuadraticShapeNCEffect.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
QuadraticShapeNCEffect::QuadraticShapeNCEffect(const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double QuadraticShapeNCEffect::calculateChangeContribution(int actor,
		int difference)
{
	return (2 * this->value(actor) + difference) * difference;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double QuadraticShapeNCEffect::egoStatistic(int ego, double * currentValues)
{
	double currentValues_nc = currentValues[ego]+this->overallCenterMean();
	return currentValues_nc * currentValues_nc;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double QuadraticShapeNCEffect::endowmentStatistic(const int * difference,
	double * currentValues)
{
	double statistic = 0;
	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
		{
			double currentValues_nc = currentValues[i]+this->overallCenterMean();
			statistic += currentValues_nc * currentValues_nc
				- (currentValues_nc + difference[i] + currentValues_nc) *
				(currentValues_nc + difference[i] + currentValues_nc);
		}
	}

	return statistic;
}

}
