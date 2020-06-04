/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreerEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeEffect class.
 *****************************************************************************/

#ifndef OUTDEGREEEFFECT_H_
#define OUTDEGREEEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Outdegree effect defined as the product of the ego with ths number of
 * its outward neighbors (with respect to a certain network).
 */
class OutdegreeEffect : public NetworkDependentBehaviorEffect
{
public:
	OutdegreeEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*OUTDEGREEEFFECT_H_*/

