/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateEgoEffect class.
 *****************************************************************************/

#ifndef COVARIATEEGOEFFECT_H_
#define COVARIATEEGOEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-ego or covariate-related activity effect (see manual).
 */
class CovariateEgoEffect : public CovariateDependentNetworkEffect
{
public:
	explicit CovariateEgoEffect(const EffectInfo * pEffectInfo,
					const bool leftThresholded, const bool rightThresholded);
	CovariateEgoEffect(const EffectInfo * pEffectInfo,
					const bool leftThresholded, const bool rightThresholded,
					const bool simulatedState);

	virtual double calculateContribution(int alter) const;
	virtual bool egoEffect() const;

protected:
	virtual double tieStatistic(int alter);


private:
	bool lleftThresholded {};
	bool lrightThresholded {};
	double lthreshold {};
};

}

#endif /*COVARIATEEGOEFFECT_H_*/
