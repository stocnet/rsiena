/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ThresholdShapeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ThresholdShapeEffect class.
 *****************************************************************************/

#include "ThresholdShapeEffect.h"
#include "model/variables/BehaviorVariable.h"
#include "model/EffectInfo.h"

namespace siena
{

/**
 * Constructor.
 */
ThresholdShapeEffect::ThresholdShapeEffect(const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo)
{
	this->lpar = pEffectInfo->internalEffectParameter();
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double ThresholdShapeEffect::calculateChangeContribution(int actor,
	int difference)
{
	double statistic = 0;
	if (((this->centeredValue(actor) + difference) >= this->lpar)
			&& (this->centeredValue(actor) < this->lpar))
	{
		statistic = 1;
	}
	else if (((this->centeredValue(actor) + difference) < this->lpar)
			&& (this->centeredValue(actor) >= this->lpar))
	{
		statistic = -1;
	}
	return statistic;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double ThresholdShapeEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	if (currentValues[ego] >= this->lpar)
	{
		statistic = 1;
	}
	return statistic;
}


/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double ThresholdShapeEffect::egoEndowmentStatistic(int ego,
	const int * difference, double * currentValues)
{
	double statistic = 0;
	if ((currentValues[ego] >= this-> lpar) && (difference[ego] > 0))
	{
		statistic = 1;
	}
	return statistic;
}

}
