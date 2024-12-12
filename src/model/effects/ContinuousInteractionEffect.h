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

#ifndef CONTINUOUSINTERACTIONEFFECT_H_
#define CONTINUOUSINTERACTIONEFFECT_H_

#include "ContinuousEffect.h"

namespace siena
{

/**
 * A user-defined interaction effect between two or three other effects.
 */
class ContinuousInteractionEffect: public ContinuousEffect
{
public:
    ContinuousInteractionEffect(const EffectInfo * pEffectInfo,
		ContinuousEffect * pEffect1,
        ContinuousEffect * pEffect2,
		ContinuousEffect * pEffect3);
	virtual ~ContinuousInteractionEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double calculateChangeContribution(int actor);

protected:
	virtual double egoStatistic(int ego,
		double * currentValues);

private:
	// The interacting effects (lpEffect3 = 0 for two-way interactions)

	ContinuousEffect * lpEffect1;
	ContinuousEffect * lpEffect2;
	ContinuousEffect * lpEffect3;
};

}

#endif /* CONTINUOUSINTERACTIONEFFECT_H_ */
