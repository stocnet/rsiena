/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * NetworkDependentBehaviorEffect class.
 *****************************************************************************/

#ifndef NETWORKDEPENDENTBEHAVIOREFFECT_H_
#define NETWORKDEPENDENTBEHAVIOREFFECT_H_

#include "BehaviorEffect.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;


// ----------------------------------------------------------------------------
// Section: NetworkDependentBehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects depending on some network variable.
 */
class NetworkDependentBehaviorEffect : public BehaviorEffect
{
public:
	NetworkDependentBehaviorEffect(const EffectInfo * pEffectInfo);
	// for gmom:
	NetworkDependentBehaviorEffect(const EffectInfo * pEffectInfo, const bool simulatedState);
	virtual ~NetworkDependentBehaviorEffect();

	virtual void initialize(const Data * pData,
			State * pState, int period, Cache * pCache);
	virtual void initialize(const Data *pData,
			State *pState, State *pSimulatedState, int period, Cache *pCache);

protected:
	inline const Network * pNetwork() const;
	double totalAlterValue(int i) const;
	double totalInAlterValue(int i) const;
	int numberAlterHigher(int i) const;
	int numberAlterLower(int i) const;
	int numberAlterEqual(int i) const;
	int numberAlterHigherPop(int i) const;
	int numberAlterLowerPop(int i) const;
	int numberAlterEqualPop(int i) const;
	virtual void preprocessEgo(int ego);

private:
	//! If `1` value(), missing() and similarity() returns the simulated value
	//! (if the covariate is a behavior) or the observed value at the end of the
	//! period.
	const int lSimulatedOffset {};

	// The network this effect is interacting with
	const Network * lpNetwork;
	// total out- and in-alter values
	double * ltotalAlterValues {};
	double * ltotalInAlterValues {};
	// number of higher, lower, and equal alter values
	int * lnumberAlterHigher {};
	int * lnumberAlterLower {};
	int * lnumberAlterEqual {};
	// and weighted by alter indegrees
	int * lnumberAlterHigherPop {};
	int * lnumberAlterLowerPop {};
	int * lnumberAlterEqualPop {};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network this effect is interacting with.
 */
const Network * NetworkDependentBehaviorEffect::pNetwork()
	const
{
	return this->lpNetwork;
}

}

#endif /*NETWORKDEPENDENTBEHAVIOREFFECT_H_*/
