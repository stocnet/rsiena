/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedNetworkEffect.h
 *
 * Description: This file contains the definition of the
 * MixedNetworkEffect class.
 *****************************************************************************/

#ifndef MIXEDNETWORKEFFECT_H_
#define MIXEDNETWORKEFFECT_H_

#include <string>
#include "NetworkEffect.h"
#include "utils/NamedObject.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class Cache;
class TwoNetworkCache;


// ----------------------------------------------------------------------------
// Section: MixedNetworkEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all network effects.
 */
class MixedNetworkEffect : public NetworkEffect
{

public:
	MixedNetworkEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	inline const Network * pNetwork() const;
	inline const NetworkLongitudinalData * pData() const;

protected:
//	bool firstInTieExists(int alter) const; // requires more in TwoNetworkCache
	bool firstOutTieExists(int alter) const;
//	bool secondInTieExists(int alter) const;
	bool secondOutTieExists(int alter) const;
	inline const Network * pFirstNetwork() const;
	inline const Network * pSecondNetwork() const;
	inline TwoNetworkCache * pTwoNetworkCache() const;

private:
	const Network * lpFirstNetwork;
	const Network * lpSecondNetwork;
	TwoNetworkCache * lpTwoNetworkCache;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

const Network * MixedNetworkEffect::pFirstNetwork() const
{
	return this->lpFirstNetwork;
}

const Network * MixedNetworkEffect::pSecondNetwork() const
{
	return this->lpSecondNetwork;
}

TwoNetworkCache * MixedNetworkEffect::pTwoNetworkCache() const
{
	return this->lpTwoNetworkCache;
}

}

#endif /*MIXEDNETWORKEFFECT_H_*/
