/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateDependentNetworkEffect.h
 *
 * Description: This file contains the definition of the
 * DyadicCovariateDependentNetworkEffect class.
 *****************************************************************************/

#ifndef DYADICCOVARIATEDEPENDENTNETWORKEFFECT_H_
#define DYADICCOVARIATEDEPENDENTNETWORKEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantDyadicCovariate;
class ChangingDyadicCovariate;
class DyadicCovariateValueIterator;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for network effects depending on a dyadic covariate.
 */
class DyadicCovariateDependentNetworkEffect : public NetworkEffect
{
public:
	DyadicCovariateDependentNetworkEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	double value(int i, int j) const;
	bool missing(int i, int j) const;
	DyadicCovariateValueIterator rowValues(int i) const;
	DyadicCovariateValueIterator columnValues(int j) const;
	bool constantCovariate() const;
	virtual void initializeStatisticCalculation();
	virtual void cleanupStatisticCalculation();

private:
	// The constant covariate this effect depends on or 0, if the
	// effect depends on a changing covariate.

	ConstantDyadicCovariate * lpConstantCovariate;

	// The changing covariate this effect depends on or 0, if the
	// effect depends on a constant covariate.

	ChangingDyadicCovariate * lpChangingCovariate;

	// flag to control exclusion of missing values
	
	bool lexcludeMissings {};
};

}

#endif /*DYADICCOVARIATEDEPENDENTNETWORKEFFECT_H_*/
