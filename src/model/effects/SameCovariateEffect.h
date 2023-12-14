/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateEffect.h
 *
 * Description: This file contains the definition of the
 * SameCovariateEffect class.
 *****************************************************************************/

#ifndef SAMECOVARIATEEFFECT_H_
#define SAMECOVARIATEEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Same covariate and same covariate x reciprocity effects (see manual).
 */
class SameCovariateEffect : public CovariateDependentNetworkEffect
{
public:
	SameCovariateEffect(const EffectInfo * pEffectInfo, 
					bool same, bool reciprocal);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the reciprocal version of the same covariate effect is
	// required

	bool lsame {};
	bool lreciprocal {};
};

}

#endif /*SAMECOVARIATEEFFECT_H_*/
