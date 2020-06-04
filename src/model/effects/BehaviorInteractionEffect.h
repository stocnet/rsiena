/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorInteractionEffect.h
 *
 * Description: This file contains the definition of the
 * BehaviorInteractionEffect class.
 *****************************************************************************/

#ifndef BEHAVIORINTERACTIONEFFECT_H_
#define BEHAVIORINTERACTIONEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * A user-defined interaction effect between two or three other effects.
 */
class BehaviorInteractionEffect: public BehaviorEffect
{
public:
	BehaviorInteractionEffect(const EffectInfo * pEffectInfo,
		BehaviorEffect * pEffect1,
		BehaviorEffect * pEffect2,
		BehaviorEffect * pEffect3);
	virtual ~BehaviorInteractionEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double calculateChangeContribution(int actor, int difference);

protected:
	virtual double egoStatistic(int ego,
		double * currentValues);
	virtual double egoEndowmentStatistic(int i, const int * difference,
		double * currentValues);

private:
	// The interacting effects (lpEffect3 = 0 for two-way interactions)

	BehaviorEffect * lpEffect1;
	BehaviorEffect * lpEffect2;
	BehaviorEffect * lpEffect3;
};

}

#endif /* BEHAVIORINTERACTIONEFFECT_H_ */
