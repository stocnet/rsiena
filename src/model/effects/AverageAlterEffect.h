/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterEffect.h
 *
 * Description: This file contains the definition of the
 * AverageAlterEffect class.
 *****************************************************************************/

#ifndef AVERAGEALTEREFFECT_H_
#define AVERAGEALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average alter effect defined as the product of the ego with the average
 * of its neighbors (with respect to a certain network).
 */
class AverageAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	AverageAlterEffect(const EffectInfo * pEffectInfo, bool divide,
		bool alterPopularity);
	AverageAlterEffect(const EffectInfo * pEffectInfo, bool divide,
		bool alterPopularity, const bool simulatedState);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	// divide indicates whether there will be division by the outdegree
	bool ldivide {};
	// alterPopularity indicates weighting by alters' indegrees
	bool lalterPopularity {};
};

}

#endif /*AVERAGEALTEREFFECT_H_*/
