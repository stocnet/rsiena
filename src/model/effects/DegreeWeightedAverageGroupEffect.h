/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DegreeWeightedAverageGroupEffect.h
 *
 * Description: This file contains the definition of the
 * DegreeWeightedAverageGroupEffect class.
 *****************************************************************************/

#ifndef DEGREEWEIGHTEDAVERAGEGROUPEFFECT_H_
#define DEGREEWEIGHTEDAVERAGEGROUPEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average of the statistic z_j weighted by in or outdegree for the group.
 */
class DegreeWeightedAverageGroupEffect : public NetworkDependentBehaviorEffect
{
public:
	DegreeWeightedAverageGroupEffect(const EffectInfo * pEffectInfo, bool divide,
	bool outdegree, bool nc);

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
	bool loutdegree {}; // true = weight by outdegree rather than indegree
	bool lnc {};
};

}

#endif /* DEGREEWEIGHTEDAVERAGEGROUPEFFECT_H_ */
