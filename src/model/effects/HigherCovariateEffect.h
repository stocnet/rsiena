/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HigherCovariateEffect.h
 *
 * Description: This file contains the definition of the
 * HigherCovariateEffect class.
 *****************************************************************************/

#ifndef HIGHERCOVARIATEEFFECT_H_
#define HIGHERCOVARIATEEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * A covariate-dependent network effect defined as
 * s_i = sum_j x_{ij} I{v_i > v_j} + 0.5 sum_j x_{ij} I{v_i = v_j}.
 */
class HigherCovariateEffect : public CovariateDependentNetworkEffect
{
public:
	HigherCovariateEffect(const EffectInfo * pEffectInfo, bool center);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool lcenter {};
};

}

#endif /*HIGHERCOVARIATEEFFECT_H_*/
