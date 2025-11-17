/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: LinearShapeNCEffect.h
 *
 * Description: This file contains the definition of the
 * LinearShapeEffect class.
 *****************************************************************************/

#ifndef LINEARSHAPENCEFFECT_H_
#define LINEARSHAPENCEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * Linear shape effect having the statistic z_i.
 */
class LinearShapeNCEffect : public BehaviorEffect
{
public:
	LinearShapeNCEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*LINEARSHAPENCEFFECT_H_*/
