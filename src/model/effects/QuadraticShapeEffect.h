/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: QuadraticShapeEffect.h
 *
 * Description: This file contains the definition of the
 * QuadraticShapeEffect class.
 *****************************************************************************/

#ifndef QUADRATICSHAPEEFFECT_H_
#define QUADRATICSHAPEEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * Quadratic shape effect having the statistic z_i^2.
 */
class QuadraticShapeEffect : public BehaviorEffect
{
public:
	QuadraticShapeEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double endowmentStatistic(const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*QUADRATICSHAPEEFFECT_H_*/
