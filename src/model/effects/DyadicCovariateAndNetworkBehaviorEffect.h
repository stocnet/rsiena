/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Behavior Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateAndNetworkBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * DyadicCovariateAndNetworkBehaviorEffect class.
 * komt van DyadicCovariateDependentNetworkEffect.h; netwerk veranderd naar behavior;
 * maar ook met invloeden DyadicCovariateDependentNetworkEffect
 *****************************************************************************/

#ifndef DYADICCOVARIATEANDNETWORKBEHAVIOREFFECT_H_
#define DYADICCOVARIATEANDNETWORKBEHAVIOREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

//class Network;
//class BehaviorVariable;
class BehaviorLongitudinalData;
class ConstantDyadicCovariate;
class ChangingDyadicCovariate;
class DyadicCovariateValueIterator;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for network effects depending on a dyadic covariate.
 */
class DyadicCovariateAndNetworkBehaviorEffect : public NetworkDependentBehaviorEffect
{
public:
	DyadicCovariateAndNetworkBehaviorEffect(const EffectInfo * pEffectInfo);
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	double dycoValue(int i, int j) const;
	bool missingDyCo(int i, int j) const;
	DyadicCovariateValueIterator rowValues(int i) const;
	DyadicCovariateValueIterator columnValues(int j) const;
	bool constantDyadicCovariate() const;
virtual void initializeStatisticCalculation();
virtual void cleanupStatisticCalculation();

private:
	// The constant covariate this effect depends on or 0, if the
	// effect depends on a changing covariate.
	ConstantDyadicCovariate * lpConstantDyadicCovariate;

	// The changing covariate this effect depends on or 0, if the
	// effect depends on a constant covariate.
	ChangingDyadicCovariate * lpChangingDyadicCovariate;

	BehaviorLongitudinalData * lpBehaviorData;
	
	// flag to control exclusion of missing values	
	bool lexcludeMissings {};
};
}

#endif /*DYADICCOVARIATEANDNETWORKBEHAVIOREFFECT_H_*/
