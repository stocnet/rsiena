/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDiffEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDiffEffect class.
 *****************************************************************************/

#ifndef COVARIATEDIFFEFFECT_H_
#define COVARIATEDIFFEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-diff effect (see manual).
 */
class CovariateDiffEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateDiffEffect(const EffectInfo * pEffectInfo, bool diff, int trafo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool ldiff {};
	bool lsquared {};
	bool labs {};
};

}

#endif /*COVARIATEDIFFEFFECT_H_*/
