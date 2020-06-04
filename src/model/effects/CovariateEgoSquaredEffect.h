/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoSquaredEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateEgoSquaredEffect class.
 *****************************************************************************/

#ifndef COVARIATEEGOSQUAREDEFFECT_H_
#define COVARIATEEGOSQUAREDEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-alter and covariate squared-alter effects (see manual).
 */
class CovariateEgoSquaredEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateEgoSquaredEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;
	virtual bool egoEffect() const;

protected:
	virtual double tieStatistic(int alter);

};

}

#endif /*COVARIATEEGOSQUAREDEFFECT_H_*/
