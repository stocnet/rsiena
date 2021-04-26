/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantEffect.h
 *
 * Description: This file contains the definition of the
 * ConstantEffect class.
 *****************************************************************************/

#ifndef CONSTANTEFFECT_H_
#define CONSTANTEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * Linear shape effect having the statistic z_i.
 */
class ConstantEffect : public BehaviorEffect
{
public:
	ConstantEffect(const EffectInfo * pEffectInfo);

	void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
	
};

}

#endif /*CONSTANTEFFECT_H_*/
