/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedOutStarFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * MixedOutStarFunction.
 *****************************************************************************/

#include "MixedOutStarFunction.h"
#include "model/tables/TwoNetworkCache.h"
#include "model/tables/MixedEgocentricConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
MixedOutStarFunction::MixedOutStarFunction(string firstNetworkName,
	string secondNetworkName) :
	MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
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
void MixedOutStarFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpTable = this->pTwoNetworkCache()->pOutStarTable();
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double MixedOutStarFunction::value(int alter)
{
	return this->lpTable->get(alter);
}


/**
 * Returns the value of this function as an integer.
 */
int MixedOutStarFunction::intValue(int alter)
{
	return this->lpTable->get(alter);
}

}
