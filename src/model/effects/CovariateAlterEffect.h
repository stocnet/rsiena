/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: https://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateAlterEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateAlterEffect class.
 *****************************************************************************/

#ifndef COVARIATEALTEREFFECT_H_
#define COVARIATEALTEREFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-alter and covariate squared-alter effects (see manual).
 */
class CovariateAlterEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateAlterEffect(const EffectInfo * pEffectInfo, const bool leftThresholded,
							const bool rightThresholded, const bool squared);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool lleftThresholded {};
	bool lrightThresholded {};
	double lthreshold {};
	bool lsquared {};
};

}

#endif /*COVARIATEALTEREFFECT_H_*/
