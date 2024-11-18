/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PopularityAlterEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeEffect class.
 *****************************************************************************/

#ifndef POPULARITYALTEREFFECT_H_
#define POPULARITYALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average and total Popularity alter effect (see manual). Also the effect of the (average) number of 
 * two-in-stars originating from ego.
 */
class PopularityAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	PopularityAlterEffect(const EffectInfo * pEffectInfo, bool divide);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	double averageInDegree(int i) const;
	// divide indicates whether there will be division by the outdegree
	bool ldivide {};
};

}

#endif /*POPULARITYALTEREFFECT_H_*/

