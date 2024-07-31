/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedNetworkAlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * MixedNetworkAlterFunction.
 *****************************************************************************/
 
#include "MixedNetworkAlterFunction.h"
#include "network/Network.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/TwoNetworkCache.h"
#include "network/CommonNeighborIterator.h"

using namespace std;

namespace siena
{

MixedNetworkAlterFunction::MixedNetworkAlterFunction(string firstNetworkName,
	string secondNetworkName )
{
	this->lfirstNetworkName = firstNetworkName;
	this->lsecondNetworkName = secondNetworkName;
	this->lpFirstNetwork = 0;
	this->lpSecondNetwork = 0;
	this->lpTwoNetworkCache = 0;
	this->lpFirstNetworkCache = 0;
}

MixedNetworkAlterFunction::~MixedNetworkAlterFunction()
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void MixedNetworkAlterFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpFirstNetwork = pState->pNetwork(this->lfirstNetworkName);
	this->lpSecondNetwork = pState->pNetwork(this->lsecondNetworkName);
	this->lpTwoNetworkCache = pCache->pTwoNetworkCache(this->lpFirstNetwork,
		this->lpSecondNetwork);
	this->lpFirstNetworkCache = pCache->pNetworkCache(this->lpFirstNetwork);
}

// ----------------------------------------------------------------------------
// Section: Tie variables
// ----------------------------------------------------------------------------

/**
 * Returns if in the first network
 * there is a tie from the current ego to the given alter.
 */
bool MixedNetworkAlterFunction::firstOutTieExists(int alter) const
{
	return this->lpTwoNetworkCache->firstOutTieExists(alter);
}

/**
 * Returns if in the second network
 * there is a tie from the current ego to the given alter.
 */
bool MixedNetworkAlterFunction::secondOutTieExists(int alter) const
{
	return this->lpTwoNetworkCache->secondOutTieExists(alter);
}

// ----------------------------------------------------------------------------
// Section: Iterators
// ----------------------------------------------------------------------------

/**
 * Returns an iterator over first-network instars of actors i and j.
 */
CommonNeighborIterator MixedNetworkAlterFunction::firstNetworkInStars(int i, int j) const 
{
//	this->lpFirstNetwork->checkSenderRange(i);
//	this->lpFirstNetwork->checkSenderRange(j);
	return CommonNeighborIterator(this->lpFirstNetwork->outTies(i), this->lpFirstNetwork->outTies(j));
}



}
