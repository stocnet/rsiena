/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageGroupEffect.h
 *
 * Description: This file contains the definition of the
 * AverageGroupEffect class.
 *****************************************************************************/

#ifndef AVERAGEGROUPEFFECT_H_
#define AVERAGEGROUPEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * Average of the statistic z_j for the group.
 */
class AverageGroupEffect : public BehaviorEffect
{
public:
	AverageGroupEffect(const EffectInfo * pEffectInfo);

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
};

}

#endif /*AVERAGEGROUPEFFECT_H_*/
