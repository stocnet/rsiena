/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageInAlterEffect.h
 *
 * Description: This file contains the definition of the
 * AverageInAlterEffect class.
 *****************************************************************************/

#ifndef AVERAGEINALTEREFFECT_H_
#define AVERAGEINALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average in-alter effect defined as the product of the ego with the average
 * of its incoming neighbors (with respect to a certain network).
 */
class AverageInAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	AverageInAlterEffect(const EffectInfo * pEffectInfo, bool divide);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	// divide indicates whether there will be division by the indegree
	bool ldivide {};
};

}

#endif /*AVERAGEINALTEREFFECT_H_*/
