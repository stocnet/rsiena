/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleDegreeBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeEffect class.
 *****************************************************************************/

#ifndef DOUBLEDEGREEBEHAVIOREFFECT_H_
#define DOUBLEDEGREEBEHAVIOREFFECT_H_

#include "TwoNetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Double (multiplex) degree behavior effects (see manual).
 */
class DoubleDegreeBehaviorEffect : public TwoNetworkDependentBehaviorEffect
{
public:
	DoubleDegreeBehaviorEffect(const EffectInfo * pEffectInfo,
				bool firstDirection, int secondDirection);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	int calculateDoubleDegree(int actor) const;
	// secondDirection can be 0, 1, 2,
	// indicating the direction of the second tie
	bool lfirstDirection {};
	int lsecondDirection {};
	// Must the degree be subtracted?
	bool lsubtract {};
};

}

#endif /*DOUBLEDEGREEBEHAVIOREFFECT_H_*/

