/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreeWeightedAverageGroupEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeWeightedAverageGroupEffect class.
 *****************************************************************************/

#ifndef INDEGREEWEIGHTEDAVERAGEGROUPEFFECT_H_
#define INDEGREEWEIGHTEDAVERAGEGROUPEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average of the statistic z_j weighted by indegree for the group.
 */
class IndegreeWeightedAverageGroupEffect : public NetworkDependentBehaviorEffect
{
public:
	IndegreeWeightedAverageGroupEffect(const EffectInfo * pEffectInfo, bool divide);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);

private:
	// lcentermean = whether to center about the general mean
	bool lcenterMean {};
	// if not lcenter, centering is about the following value
	double lcenteringValue {};
	bool ldivide {};
};

}

#endif /*INDEGREEWEIGHTEDAVERAGEGROUPEFFECT_H_*/
