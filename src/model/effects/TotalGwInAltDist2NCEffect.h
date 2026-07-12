/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalGwInAltDist2NCEffect.h
 *
 * Description: This file contains the definition of the
 * TotalGwInAltDist2NCEffect class.
 *****************************************************************************/

#ifndef TOTALGWINALTDIST2NCEFFECT_H_
#define TOTALGWINALTDIST2NCEFFECT_H_

#include <vector>
#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Non-centered alter-wise geometrically weighted in-distance-2 effect.
 */
class TotalGwInAltDist2NCEffect : public NetworkDependentBehaviorEffect
{
public:
	explicit TotalGwInAltDist2NCEffect(const EffectInfo * pEffectInfo, bool forward);

	virtual void initialize(const Data * pData, State * pState, int period,
		Cache * pCache) override;
	virtual double calculateChangeContribution(int actor,
		int difference) override;
	virtual double egoStatistic(int ego, double * currentValues) override;

private:
	std::vector<double> lcumulativeWeight;

	double lparameter {};
	double lweight {};
	double lexpmweight {};
	double lexpfactor {};
	bool lforward {};
};

}

#endif /* TOTALGWINALTDIST2NCEFFECT_H_ */