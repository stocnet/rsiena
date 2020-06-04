/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IsolateOutContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeEffect class.
 *****************************************************************************/

#ifndef ISOLATEOUTCONTINUOUSEFFECT_H_
#define ISOLATEOUTCONTINUOUSEFFECT_H_

#include "NetworkDependentContinuousEffect.h"

namespace siena
{

/**
 * Isolate oudegree effect defined as 1 in case an actor has no outgoing ties 
 * and 0 otherwise (with respect to a certain network).
 */
class IsolateOutContinuousEffect : public NetworkDependentContinuousEffect
{
public:
	IsolateOutContinuousEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*ISOLATEOUTCONTINUOUSEFFECT_H_*/

