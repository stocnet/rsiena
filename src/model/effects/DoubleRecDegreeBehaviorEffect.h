/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleRecDegreeBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeEffect class.
 *****************************************************************************/

#ifndef DOUBLERECDEGREEBEHAVIOREFFECT_H_
#define DOUBLERECDEGREEBEHAVIOREFFECT_H_

#include "TwoNetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Double (multiplex) degree behavior effects (see manual).
 */
class DoubleRecDegreeBehaviorEffect : public TwoNetworkDependentBehaviorEffect
{
public:
	DoubleRecDegreeBehaviorEffect(const EffectInfo * pEffectInfo,
				int secondDirection);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	int calculateDoubleRecDegree(int actor) const;
	// secondDirection indicates the direction of the second tie
	int lsecondDirection {};
};

}

#endif /*DOUBLERECDEGREEBEHAVIOREFFECT_H_*/

