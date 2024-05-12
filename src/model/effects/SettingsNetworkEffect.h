/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SettingsNetworkEffect.h
 *
 * Description: This file contains the definition of the
 * SettingsNetworkEffect class.
 *****************************************************************************/

#ifndef SETTINGSNETWORKEFFECT_H_
#define SETTINGSNETWORKEFFECT_H_

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
// Section: SettingsNetworkEffect class
// ----------------------------------------------------------------------------

/**
 * Network effects depending on primary setting
 */
class SettingsNetworkEffect : public NetworkEffect
{

public:
	SettingsNetworkEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	inline const Network * pNetwork() const;
	virtual void preprocessEgo(int ego);

protected:
	bool outTieExists(int alter) const;
	bool settingsTieExists(int alter) const;
	int outDegree() const;
	int settingDegree() const;
	inline const Network * pPrimarySetting() const;
	inline TwoNetworkCache * pTwoNetworkCache() const;
	inline int stepType() const;

private:
	const Network * lpNetwork;
	const Network * lpPrimarySetting;
	TwoNetworkCache * lpTwoNetworkCache;
	int lstepType {};

};
// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

const Network * SettingsNetworkEffect::pNetwork() const
{
	return this->lpNetwork;
}

int SettingsNetworkEffect::stepType() const
{
	return this->lstepType;
}

const Network * SettingsNetworkEffect::pPrimarySetting() const
{
	return this->lpPrimarySetting;
}

TwoNetworkCache * SettingsNetworkEffect::pTwoNetworkCache() const
{
	return this->lpTwoNetworkCache;
}

}

#endif /*SETTINGSNETWORKEFFECT_H_*/
