/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocalDegreeBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeEffect class.
 *****************************************************************************/

#ifndef RECIPROCALDEGREEBEHAVIOREFFECT_H_
#define RECIPROCALDEGREEBEHAVIOREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Reciprocated degree behavior effect (see manual).
 */
class ReciprocalDegreeBehaviorEffect : public NetworkDependentBehaviorEffect
{
public:
	ReciprocalDegreeBehaviorEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*RECIPROCALDEGREEBEHAVIOREFFECT_H_*/

