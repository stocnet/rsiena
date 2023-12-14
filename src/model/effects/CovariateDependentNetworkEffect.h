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

#ifndef COVARIATEDEPENDENTNETWORKEFFECT_H_
#define COVARIATEDEPENDENTNETWORKEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;
class ChangingCovariate;
class BehaviorVariable;
class BehaviorLongitudinalData;
class ContinuousLongitudinalData;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for network effects depending on an individual
 * covariate (constant, changing, or dependent behavior variable).
 */
class CovariateDependentNetworkEffect : public NetworkEffect
{
public:
	explicit CovariateDependentNetworkEffect(const EffectInfo * pEffectInfo);
	CovariateDependentNetworkEffect(const EffectInfo * pEffectInfo, const bool simulated);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void initialize(const Data * pData, State * pState,
			State * pSimulatedState, int period, Cache * pCache);

protected:
	double value(const int i) const;
	bool missing(int i) const;
	double actor_similarity(int i, int j) const;
	ConstantCovariate * pConstantCovariate() const;
	ChangingCovariate * pChangingCovariate() const;
	BehaviorLongitudinalData * pBehaviorData() const;
    ContinuousLongitudinalData * pContinuousData() const;

private:
	//! If `1` value(), missing() and actor_similarity() returns the simulated value
	//! (if the covariate is a behavior) or the observed value at the end of the
	//! period.
	const int lSimulatedOffset {};

	ConstantCovariate * lpConstantCovariate;
	ChangingCovariate * lpChangingCovariate;
	BehaviorLongitudinalData * lpBehaviorData;
    ContinuousLongitudinalData * lpContinuousData;

	// The current value of a (discrete or continuous) behavior
    // variable per each actor.
	// This array is 0 for covariate-based effects.

	const int * lvalues {};
   const double * lcontinuousValues {};
};

}

#endif /*COVARIATEDEPENDENTNETWORKEFFECT_H_*/
