/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * MixedNetworkEffect class.
 *****************************************************************************/

#include "MixedNetworkEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/Cache.h"
#include "model/tables/TwoNetworkCache.h"
#include "model/EffectInfo.h"
#include "model/State.h"



using namespace std;

namespace siena
{

/**
 * Constructor.
 */
MixedNetworkEffect::MixedNetworkEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
	this->lpFirstNetwork = 0;
	this->lpSecondNetwork = 0;
	this->lpTwoNetworkCache = 0;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void MixedNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	string networkName1 = this->pEffectInfo()->variableName();
	string networkName2 = this->pEffectInfo()->interactionName1();
	this->lpFirstNetwork = pState->pNetwork(networkName1);
	this->lpSecondNetwork = pState->pNetwork(networkName2);
	this->lpTwoNetworkCache = pCache->pTwoNetworkCache(this->lpFirstNetwork,
		this->lpSecondNetwork);
}

/**
 * Returns if in the first network
 * there is a tie from the current ego to the given alter.
 */
bool MixedNetworkEffect::firstOutTieExists(int alter) const
{
	return this->lpTwoNetworkCache->firstOutTieExists(alter);
}

/**
 * Returns if in the second network
 * there is a tie from the current ego to the given alter.
 */
bool MixedNetworkEffect::secondOutTieExists(int alter) const
{
	return this->lpTwoNetworkCache->secondOutTieExists(alter);
}


}
