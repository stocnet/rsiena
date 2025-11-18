/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: QuadraticShapeCcEffect.cpp
 *
 * Description: This file contains the implementation of the
 * QuadraticShapeCcEffect class.
 *****************************************************************************/

#include "QuadraticShapeCcEffect.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
QuadraticShapeCcEffect::QuadraticShapeCcEffect(const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double QuadraticShapeCcEffect::calculateChangeContribution(int actor,
		int difference)
{
	double currentMean = 0; // contemporaneous mean, inefficiently calculated each time here

	for (int i = 0; i < this->n(); i++)
	{
		currentMean += this->value(i);
	}
	currentMean /= this->n();

	return (2 * this->value(actor) - currentMean + difference) * difference;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double QuadraticShapeCcEffect::egoStatistic(int ego, double * currentValues)
{
	double currentDelta = 0; // correction for contemporaneous mean, inefficiently calculated each time here

	for (int j = 0; j < this->n(); j++)
	{
		currentDelta += currentValues[j]; // note these are centered at grand mean
	}
	currentDelta /= this->n();

	return (currentValues[ego]-currentDelta) * (currentValues[ego] - currentDelta);
}


}
