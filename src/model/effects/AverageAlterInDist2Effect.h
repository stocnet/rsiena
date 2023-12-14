/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterInDist2Effect.h
 *
 * Description: This file contains the definition of the
 * AverageAlterInDist2Effect class.
 *****************************************************************************/

#ifndef AVERAGEALTERINDIST2EFFECT_H_
#define AVERAGEALTERINDIST2EFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average alter effect defined as the product of the ego with the average
 * of the alter average of its neighbors (with respect to a certain network).
 */
class AverageAlterInDist2Effect : public NetworkDependentBehaviorEffect
{
public:
	AverageAlterInDist2Effect(const EffectInfo * pEffectInfo,
						bool divide1, bool divide2);
	virtual ~AverageAlterInDist2Effect();
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	bool ldivide1 {};
	// Indicates whether there will be division by the outdegree of ego
	bool ldivide2 {};
	// Indicates whether there will be division by the indegree of alter
};

}

#endif /*AVERAGEALTERINDIST2EFFECT_H_*/
