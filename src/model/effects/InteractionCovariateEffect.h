/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InteractionCovariateEffect.h
 *
 * Description: This file contains the definition of the
 * InteractionCovariateEffect class.
 *****************************************************************************/

#ifndef INTERACTIONCOVARIATEEFFECT_H_
#define INTERACTIONCOVARIATEEFFECT_H_

#include "CovariateDependentBehaviorEffect.h"

namespace siena
{

/**
 * Interaction covariate behavior effect (see manual).
 */
class InteractionCovariateEffect : public CovariateDependentBehaviorEffect
{
public:
	InteractionCovariateEffect(const EffectInfo * pEffectInfo,
		bool averageSimilarity,
		bool totalSimilarity,
		bool averageAlter,
		bool totalAlter);
	virtual ~InteractionCovariateEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);

private:
	EffectInfo * lpInteractingEffectInfo;
	BehaviorEffect * lpInteractingEffect;
};

}

#endif /*INTERACTIONCOVARIATEEFFECT_H_*/
