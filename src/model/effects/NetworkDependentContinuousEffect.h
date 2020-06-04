/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * NetworkDependentContinuousEffect class.
 *****************************************************************************/

#ifndef NETWORKDEPENDENTCONTINUOUSEFFECT_H_
#define NETWORKDEPENDENTCONTINUOUSEFFECT_H_

#include "ContinuousEffect.h"

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
 * The base class for all continuous behavior effects depending on some network 
 * variable.
 */
class NetworkDependentContinuousEffect : public ContinuousEffect
{
public:
	NetworkDependentContinuousEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	inline const Network * pNetwork() const;

private:
	// The network this effect is interacting with
	const Network * lpNetwork;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network this effect is interacting with.
 */
const Network * NetworkDependentContinuousEffect::pNetwork()
	const
{
	return this->lpNetwork;
}

}

#endif /*NETWORKDEPENDENTCONTINUOUSEFFECT_H_*/
