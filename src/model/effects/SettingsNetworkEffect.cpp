/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SettingsNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * SettingsNetworkEffect class.
 *****************************************************************************/

#include <stdexcept>
#include "SettingsNetworkEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/TwoNetworkCache.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
SettingsNetworkEffect::SettingsNetworkEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
	this->lpNetwork = 0;
	this->lpPrimarySetting = 0;
	this->lpTwoNetworkCache = 0;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void SettingsNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	string networkName = this->pEffectInfo()->variableName();
	string pname = "primary(" + this->pEffectInfo()->variableName() + ")";
	this->lpNetwork = pState->pNetwork(networkName);
	this->lpPrimarySetting  = pState->pNetwork(pname);
	if (!lpPrimarySetting)
	{
		throw logic_error("Settings network '" +
			pname + "' expected but not found.");
	}
	this->lpTwoNetworkCache = pCache->pTwoNetworkCache(this->lpNetwork,
		this->lpPrimarySetting);
//	this->lpTwoNetworkCache->initialize(this->lego);
	this->lstepType = pCache->pNetworkCache(this->lpNetwork)->stepTypeValue();
}



/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego.
 */
void SettingsNetworkEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);
	this->lpTwoNetworkCache->initialize(ego);
//	this->lpPrimarySetting -> initSetting(ego);
	// should there be a terminateSetting somewhere?
}



/**
 * Returns if in the base network
 * there is a tie from the current ego to the given alter.
 */
bool SettingsNetworkEffect::outTieExists(int alter) const
{
	return this->lpTwoNetworkCache->firstOutTieExists(alter);
}

/**
 * Returns if in the primary setting
 * there is a tie from the current ego to the given alter.
 */
bool SettingsNetworkEffect::settingsTieExists(int alter) const
{
	return this->lpTwoNetworkCache->secondOutTieExists(alter);
}


/**
 * Returns the outdegree of the current ego
 */
int SettingsNetworkEffect::outDegree() const
{
	return this->lpNetwork->outDegree(this->ego());
}

/**
 * Returns the primary setting size (outdegree) of the current ego
 */
int SettingsNetworkEffect::settingDegree() const
{
// both of the following are OK for targets
	//return this->lpTwoNetworkCache->secondOutDegree(); 
	return this->lpPrimarySetting->outDegree(this->ego()); 
}

}
