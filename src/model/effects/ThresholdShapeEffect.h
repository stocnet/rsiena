/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ThresholdShapeEffect.h
 *
 * Description: This file contains the definition of the
 * ThresholdShapeEffect class.
 *****************************************************************************/

#ifndef THRESHOLDSHAPEEFFECT_H_
#define THRESHOLDSHAPEEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * Threshold shape effect having the statistic z_i >= par.
 */
class ThresholdShapeEffect : public BehaviorEffect
{
public:
	ThresholdShapeEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);

private:
	// lpar is the internal effect parameter
	int lpar {};
};

}

#endif /*THRESHOLDSHAPEEFFECT_H_*/
