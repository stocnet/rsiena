/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * AverageAlterContinuousEffect class.
 *****************************************************************************/

#ifndef AVERAGEALTERCONTINUOUSEFFECT_H_
#define AVERAGEALTERCONTINUOUSEFFECT_H_

#include "NetworkDependentContinuousEffect.h"

namespace siena
{

/**
 * Average alter effect defined as the average of an ego's neighbors (with
 * respect to a certain network).
 */
class AverageAlterContinuousEffect : public NetworkDependentContinuousEffect
{
public:
	AverageAlterContinuousEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*AVERAGEALTERCONTINUOUSEFFECT_H_*/
