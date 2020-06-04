/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutJaccardFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * OutJaccardFunction.
 *****************************************************************************/

#include "OutJaccardFunction.h"
#include "network/Network.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
OutJaccardFunction::OutJaccardFunction(string networkName) :
	NetworkAlterFunction(networkName)
{
	this->lpTable = 0;
}


/**
 * Outitializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void OutJaccardFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpTable = this->pNetworkCache()->pInStarTable();
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double OutJaccardFunction::value(int alter)
{
	const Network * pNetwork = this->pNetwork();
	double statistic = 0;
	int instar = this->lpTable->get(alter);
	double denom = pNetwork->outDegree(alter) +
					pNetwork->outDegree(this->ego()) - instar;
	if (denom >= 1)
	{
		statistic = instar/denom;
	}
	return statistic;
}

}
