/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WienerEffect.h
 *
 * Description: This file contains the definition of the
 * WienerEffect class.
 *****************************************************************************/

#ifndef WIENEREFFECT_H_
#define WIENEREFFECT_H_

#include "CovariateDependentContinuousEffect.h"

namespace siena
{

/**
 * Feedback effect in stochastic differential equation.
 */
class WienerEffect : public CovariateDependentContinuousEffect
{
public:
	WienerEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);

};

}

#endif /*FEEDBACKEFFECT_H_*/
