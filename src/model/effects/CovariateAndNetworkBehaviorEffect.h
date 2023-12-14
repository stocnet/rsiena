/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateAndNetworkBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateAndNetworkBehaviorEffect class.
 *****************************************************************************/

#ifndef COVARIATEANDNETWORKBEHAVIOREFFECT_H_
#define COVARIATEANDNETWORKBEHAVIOREFFECT_H_

#include "CovariateDependentBehaviorEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for behavior effects depending on an individual
 * covariate (constant, changing, or other dependent behavior variable)
 * and on a network too.
 */
class CovariateAndNetworkBehaviorEffect :
public CovariateDependentBehaviorEffect
{
public:
	CovariateAndNetworkBehaviorEffect(const EffectInfo * pEffectInfo);
	virtual ~CovariateAndNetworkBehaviorEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

protected:
	bool missingDummy(int i) const;
	bool missingInDummy(int i) const;
	double averageAlterValue(int i) const;
	double minimumAlterValue(int i) const;
	double maximumAlterValue(int i) const;
	double totalAlterValue(int i) const;
	double averageInAlterValue(int i) const;
	double totalInAlterValue(int i) const;
	inline const Network * pNetwork() const;

private:
	// The network this effect is interacting with
	const Network * lpNetwork;
	double * laverageAlterValues {};
	double * lminimumAlterValues {};
	double * lmaximumAlterValues {};
	double * ltotalAlterValues {};
	double * laverageInAlterValues {};
	double * ltotalInAlterValues {};
	bool * laverageAlterMissing {};
	bool * laverageInAlterMissing {};
};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network this effect is interacting with.
 */
const Network * CovariateAndNetworkBehaviorEffect::pNetwork()
	const
{
	return this->lpNetwork;
}
}

#endif /*COVARIATEANDNETWORKBEHAVIOREFFECT_H_*/
