/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MaxAlterContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * MaxAlterContinuousEffect class.
 *****************************************************************************/

#ifndef MAXALTERCONTINUOUSEFFECT_H_
#define MAXALTERCONTINUOUSEFFECT_H_

#include "NetworkDependentContinuousEffect.h"

namespace siena
{

/**
 * Max alter effect defined as the max of an ego's neighbors (with 
 * respect to a certain network); or min in the minim case.
 */
class MaxAlterContinuousEffect : public NetworkDependentContinuousEffect
{
public:
	MaxAlterContinuousEffect(const EffectInfo * pEffectInfo, bool minim);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	// minim indicates whether it will be a minimum instead of a maximum
	bool lminim {};
};

}

#endif /*MAXALTERCONTINUOUSEFFECT_H_*/
