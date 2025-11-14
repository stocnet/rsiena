/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterCcEffect.h
 *
 * Description: This file contains the definition of the
 * AverageAlterCcEffect class.
 *****************************************************************************/

#ifndef AVERAGEALTERCCEFFECT_H_
#define AVERAGEALTERCCEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average alter effect defined as the product of the ego with the average
 * of its neighbors (with respect to a certain network).
 */
class AverageAlterCcEffect : public NetworkDependentBehaviorEffect
{
public:
	AverageAlterCcEffect(const EffectInfo * pEffectInfo, bool divide);
	AverageAlterCcEffect(const EffectInfo * pEffectInfo, bool divide,
		const bool simulatedState);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	// divide indicates whether there will be division by the outdegree
	bool ldivide {};
};

}

#endif /*AVERAGEALTERCCEFFECT_H_*/
