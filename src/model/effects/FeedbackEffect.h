/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: FeedbackEffect.h
 *
 * Description: This file contains the definition of the
 * FeedbackEffect class.
 *****************************************************************************/

#ifndef FEEDBACKEFFECT_H_
#define FEEDBACKEFFECT_H_

#include "CovariateDependentContinuousEffect.h"

namespace siena
{

/**
 * Feedback effect in stochastic differential equation.
 */
class FeedbackEffect : public CovariateDependentContinuousEffect
{
public:
	FeedbackEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);

};

}

#endif /*FEEDBACKEFFECT_H_*/
