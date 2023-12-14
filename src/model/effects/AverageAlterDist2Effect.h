/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterDist2Effect.h
 *
 * Description: This file contains the definition of the
 * AverageAlterDist2Effect class.
 *****************************************************************************/

#ifndef AVERAGEALTERDIST2EFFECT_H_
#define AVERAGEALTERDIST2EFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average alter effect defined as the product of the ego with the average
 * of the alter average of its neighbors (with respect to a certain network).
 */
class AverageAlterDist2Effect : public NetworkDependentBehaviorEffect
{
public:
	AverageAlterDist2Effect(const EffectInfo * pEffectInfo,
						bool divide1, bool divide2);
	virtual ~AverageAlterDist2Effect();
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	bool ldivide1 {};
	// Indicates whether there will be division by the outdegree of ego
	bool ldivide2 {};
	// Indicates whether there will be division by the outdegree of alter
};

}

#endif /*AVERAGEALTERDIST2EFFECT_H_*/
