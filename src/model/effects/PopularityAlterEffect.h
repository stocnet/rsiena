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
 * Popularity alter effect (see manual).
 */
class PopularityAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	PopularityAlterEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	double averageInDegree(int i) const;
};

}

#endif /*POPULARITYALTEREFFECT_H_*/

