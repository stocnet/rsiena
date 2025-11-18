/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: QuadraticShapeCcEffect.h
 *
 * Description: This file contains the definition of the
 * QuadraticShapeCcEffect class.
 ******************	***********************************************************/

#ifndef QUADRATICSHAPECCEFFECT_H_
#define QUADRATICSHAPECCEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * Quadratic shape effect having the statistic z_i^2.
 */
class QuadraticShapeCcEffect : public BehaviorEffect
{
public:
	QuadraticShapeCcEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*QUADRATICSHAPECCEFFECT_H_*/
