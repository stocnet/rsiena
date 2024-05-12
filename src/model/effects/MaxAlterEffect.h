/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MaxAlterEffect.h
 *
 * Description: This file contains the definition of the
 * MaxAlterEffect class.
 *****************************************************************************/

#ifndef MAXALTEREFFECT_H_
#define MAXALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Max alter effect defined as the product of the ego with the max
 * of its neighbors (with respect to a certain network); or min in the minim case..
 */
class MaxAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	MaxAlterEffect(const EffectInfo * pEffectInfo, bool minim);
	MaxAlterEffect(const EffectInfo * pEffectInfo, bool minim, const bool simulatedState);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	// minim indicates whether it will be a minimum instead of a maximum
	bool lminim {};
};

}

#endif /*MAXALTEREFFECT_H_*/
