/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DegreeMixedPopularityEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeEffect class.
 *****************************************************************************/

#ifndef DEGREEMIXEDPOPULARITYEFFECT_H_
#define DEGREEMIXEDPOPULARITYEFFECT_H_

#include "TwoNetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Double (multiplex) degree behavior effects (see manual).
 */
class DegreeMixedPopularityEffect : public TwoNetworkDependentBehaviorEffect
{
public:
	DegreeMixedPopularityEffect(const EffectInfo * pEffectInfo,
				bool direction);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	int calculateMixedPopDegree(int actor) const;
	bool ldirection {};
};

}

#endif /*DEGREEMIXEDPOPULARITYEFFECT_H_*/

