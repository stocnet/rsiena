/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDependentNetworkEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDependentNetworkEffect class.
 *****************************************************************************/

#ifndef CATCOVARIATEDEPENDENTNETWORKEFFECT_H_
#define CATCOVARIATEDEPENDENTNETWORKEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for network effects depending on an individual
 * categorical covariate (constant, changing, or dependent behavior variable).
 * Covariates should be on the second mode, if network is two-mode!
 */
class CatCovariateDependentNetworkEffect : public CovariateDependentNetworkEffect
{
public:
	explicit CatCovariateDependentNetworkEffect(const EffectInfo * pEffectInfo);
	CatCovariateDependentNetworkEffect(const EffectInfo * pEffectInfo, const bool simulated);
	virtual ~CatCovariateDependentNetworkEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void initialize(const Data * pData, State * pState,
			State * pSimulatedState, int period, Cache * pCache);

protected:
	int numberCovariateTies(int a) const;
	int covariateIntValue(int i) const;

private:
	//! If `1` value(), missing() and actor_similarity() returns the simulated value
	//! (if the covariate is a behavior) or the observed value at the end of the
	//! period.
	const int lSimulatedOffset {};
	int * lpCovariateNumbers {};
};

}

#endif /*CATCOVARIATEDEPENDENTNETWORKEFFECT_H_*/
