/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageTwoInStarAlterEffect.h
 *
 * Description: This file contains the definition of the
 * AverageTwoInStarAlterEffect class.
 *****************************************************************************/

#ifndef AVERAGETWOINSTARALTEREFFECT_H_
#define AVERAGETWOINSTARALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Total two-in-star alter effect defined as the product of the behavior of ego i with the behavior
 * of the in-alters j of it's neighbors h weighted by the number of shared outstars i → h ← j 
 * (with respect to a certain network). 
 * This is a measure of influence by (locally) structural similar others.
 */
class AverageTwoInStarAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	AverageTwoInStarAlterEffect(const EffectInfo * pEffectInfo,
						bool divide1, bool divide2);
	virtual ~AverageTwoInStarAlterEffect();
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

#endif /*AVERAGETWOINSTARALTEREFFECT_H_*/
