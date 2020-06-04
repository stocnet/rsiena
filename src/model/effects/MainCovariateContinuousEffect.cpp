/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MainCovariateContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * MainCovariateContinuousEffect class.
 *****************************************************************************/

#include "MainCovariateContinuousEffect.h"

namespace siena
{

/**
 * Constructor.
 */
MainCovariateContinuousEffect::MainCovariateContinuousEffect(
	const EffectInfo * pEffectInfo) :
		CovariateDependentContinuousEffect(pEffectInfo)
{
}


/**
 * Returns how much this effect contributes to the change in the
 * continuous behavior.
 */
double MainCovariateContinuousEffect::calculateChangeContribution(int actor)
{
	return this->covariateValue(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double MainCovariateContinuousEffect::egoStatistic(int ego, double * currentValues)
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
 *//*
double MainCovariateContinuousEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if (difference[ego] > 0)
	{
		//statistic = currentValues[ego] * this->covariateValue(ego);
		statistic = - difference[ego] * this->covariateValue(ego);
	}

	return statistic;
}*/

}
