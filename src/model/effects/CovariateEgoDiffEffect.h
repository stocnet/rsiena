/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoDiffEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateEgoDiffEffect class.
 *****************************************************************************/

#ifndef COVARIATEEGODIFFEFFECT_H_
#define COVARIATEEGODIFFEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-related absolute degree difference effects (see manual).
 */
class CovariateEgoDiffEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateEgoDiffEffect(const EffectInfo * pEffectInfo,
					const bool plus, const bool minus);

	virtual double calculateContribution(int alter) const;
	virtual double endowmentStatistic(Network * pLostTieNetwork);
	virtual bool egoEffect() const;

protected:
	virtual double egoStatistic(int ego, const Network * pNetwork);

private:
	bool lplus {};
	bool lminus {};
};

}

#endif /*COVARIATEEGODIFFEFFECT_H_*/
