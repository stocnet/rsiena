/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InJaccardFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * InJaccardFunction.
 *****************************************************************************/

#include "InJaccardFunction.h"
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
InJaccardFunction::InJaccardFunction(string networkName) :
	NetworkAlterFunction(networkName)
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
void InJaccardFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpTable = this->pNetworkCache()->pOutStarTable();
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double InJaccardFunction::value(int alter)
{
	const Network * pNetwork = this->pNetwork();
	double statistic = 0;
	int outstar = this->lpTable->get(alter);
	double denom = pNetwork->inDegree(alter) +
					pNetwork->inDegree(this->ego()) - outstar;
	if (denom >= 1)
	{
		statistic = outstar/denom;
	}
	return statistic;
}

}
