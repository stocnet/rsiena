/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InterceptEffect.h
 *
 * Description: This file contains the definition of the
 * InterceptEffect class.
 *****************************************************************************/

#ifndef INTERCEPTEFFECT_H_
#define INTERCEPTEFFECT_H_

#include "ContinuousEffect.h"

namespace siena
{

/**
 * Intercept effect having the statistic ......
 */
class InterceptEffect : public ContinuousEffect
{
public:
	InterceptEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor);
//	virtual double egoEndowmentStatistic(int ego, const int * difference,
//		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

};

}

#endif /*INTERCEPTEFFECT_H_*/
