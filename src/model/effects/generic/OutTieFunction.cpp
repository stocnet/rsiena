/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutTieFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * OutTieFunction.
 *****************************************************************************/

#include <R_ext/Print.h>
#include "OutTieFunction.h"
#include "model/tables/NetworkCache.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
OutTieFunction::OutTieFunction(string networkName) :
	NetworkAlterFunction(networkName)
{
}


/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
OutTieFunction::OutTieFunction(string networkName, 
				const bool simulatedState) :
	NetworkAlterFunction(networkName, simulatedState)
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void OutTieFunction::initialize(const Data * pData,
	State * pState,  
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] pSimulatedState the current simulated state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void OutTieFunction::initialize(const Data * pData,
	State * pState, State * pSimulatedState,  
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, pSimulatedState, period, pCache);	
//	Rprintf("Passeert hier\n");
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double OutTieFunction::value(int alter) const
{
	return this->pNetworkCache()->outTieValue(alter);
}

}
