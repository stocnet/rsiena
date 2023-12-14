/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocalDegreeContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * ReciprocalDegreeContinuousEffect class.
 *****************************************************************************/

#ifndef RECIPROCALDEGREECONTINUOUSEFFECT_H_
#define RECIPROCALDEGREECONTINUOUSEFFECT_H_

#include "NetworkDependentContinuousEffect.h"

namespace siena
{

/**
 * Reciprocated degree behavior effect (see manual).
 */
class ReciprocalDegreeContinuousEffect : public NetworkDependentContinuousEffect
{
public:
	ReciprocalDegreeContinuousEffect(const EffectInfo * pEffectInfo, bool recip);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);
	
private:
	// recip indicates whether the effect addresses reciprocated or 
	// non-reciprocated ties
	bool lrecip {};
};

}

#endif /*RECIPROCALDEGREECONTINUOUSEFFECT_H_*/

