/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDiffEgoEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDiffEffect class.
 *****************************************************************************/

#ifndef COVARIATEDIFFEGOEFFECT_H_
#define COVARIATEDIFFEGOEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-diff ego effect (see manual).
 */
class CovariateDiffEgoEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateDiffEgoEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*COVARIATEDIFFEGOEFFECT_H_*/
