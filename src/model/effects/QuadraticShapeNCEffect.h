/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: QuadraticShapeNCEffect.h
 *
 * Description: This file contains the definition of the
 * QuadraticShapeNCEffect class.
 *****************************************************************************/

#ifndef QUADRATICSHAPENCEFFECT_H_
#define QUADRATICSHAPEEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * Quadratic shape effect having the statistic z_i^2.
 */
class QuadraticShapeNCEffect : public BehaviorEffect
{
public:
	QuadraticShapeNCEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double endowmentStatistic(const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*QUADRATICSHAPENCEFFECT_H_*/
