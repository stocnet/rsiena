/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MainCovariateContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * MainCovariateContinuousEffect class.
 *****************************************************************************/

#ifndef MAINCOVARIATECONTINUOUSEFFECT_H_
#define MAINCOVARIATECONTINUOUSEFFECT_H_

#include "CovariateDependentContinuousEffect.h"

namespace siena
{

/**
 * Main covariate behavior effect defined as the product of the ego
 * with the covariate.
 */
class MainCovariateContinuousEffect : public CovariateDependentContinuousEffect
{
public:
	MainCovariateContinuousEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor);
	//virtual double egoEndowmentStatistic(int ego, const int * difference,
	//	double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*MAINCOVARIATECONTINUOUSEFFECT_H_*/
