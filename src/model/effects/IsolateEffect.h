/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IsolateEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeEffect class.
 *****************************************************************************/

#ifndef ISOLATEEFFECT_H_
#define ISOLATEEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Isolate effect defined for in as z_i * I{x_{+i} = 0} and else z_i * I{x_{i+} = 0}
 */
class IsolateEffect : public NetworkDependentBehaviorEffect
{
public:
	IsolateEffect(const EffectInfo * pEffectInfo, bool in);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	// Indicates if in- or out-isolate
	bool lin {};

};

}

#endif /*ISOLATEEFFECT_H_*/

