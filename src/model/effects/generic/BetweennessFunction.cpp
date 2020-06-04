/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BetweennessFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * BetweennessFunction.
 *****************************************************************************/

#include "BetweennessFunction.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/ConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
BetweennessFunction::BetweennessFunction(string networkName) :
	OneModeNetworkAlterFunction(networkName)
{
	this->lpTable = 0;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void BetweennessFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	OneModeNetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpTable = this->pNetworkCache()->pBetweennessTable();
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before.
 */
double BetweennessFunction::value(int alter)
{
	return this->lpTable->get(alter);
}


/**
 * Returns the value of this function as an integer.
 */
int BetweennessFunction::intValue(int alter)
{
	return this->lpTable->get(alter);
}

}
