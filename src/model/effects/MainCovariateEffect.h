/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MainCovariateEffect.h
 *
 * Description: This file contains the definition of the
 * MainCovariateEffect class.
 *****************************************************************************/

#ifndef MAINCOVARIATEEFFECT_H_
#define MAINCOVARIATEEFFECT_H_

#include "CovariateDependentBehaviorEffect.h"

namespace siena
{

/**
 * Main covariate behavior effect defined as the product of the ego
 * with the covariate.
 */
class MainCovariateEffect : public CovariateDependentBehaviorEffect
{
public:
	MainCovariateEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*MAINCOVARIATEEFFECT_H_*/
