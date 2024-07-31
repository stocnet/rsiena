/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkAlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * NetworkAlterFunction.
 *****************************************************************************/

#include "NetworkAlterFunction.h"
#include "network/Network.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"

using namespace std;

namespace siena
{

NetworkAlterFunction::NetworkAlterFunction(string networkName) :
	NamedObject(networkName)
{
	this->lpNetwork = 0;
	this->lnetworkName = networkName;
	this->lpNetworkCache = 0;
}


NetworkAlterFunction::~NetworkAlterFunction()
{
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void NetworkAlterFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpNetwork = pState->pNetwork(this->lnetworkName);
	this->lpNetworkCache = pCache->pNetworkCache(this->lpNetwork);
}

/**
 * Returns if there is a tie from the current ego to the given alter.
 */
bool NetworkAlterFunction::outTieExists(int alter) const
{
	return this->lpNetworkCache->outTieExists(alter);
}

/**
 * Returns if there is a tie from the given alter to the current ego.
 */
bool NetworkAlterFunction::inTieExists(int alter) const
{
	return this->lpNetworkCache->inTieExists(alter);
}

}
