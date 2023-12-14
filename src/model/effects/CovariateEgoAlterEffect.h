/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoAlterEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateEgoAlterEffect class.
 *****************************************************************************/

#ifndef COVARIATEEGOALTEREFFECT_H_
#define COVARIATEEGOALTEREFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-ego x alter effect and covariate-ego x alter x reciprocity
 * effect (see manual).
 */
class CovariateEgoAlterEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateEgoAlterEffect(const EffectInfo * pEffectInfo, bool reciprocal);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the reciprocal version of the effect is required
	bool lreciprocal {};
};

}

#endif /*COVARIATEEGOALTEREFFECT_H_*/
