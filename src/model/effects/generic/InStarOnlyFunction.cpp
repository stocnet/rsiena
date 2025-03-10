/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InStarOnlyFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * InStarOnlyFunction.
 *****************************************************************************/

#include "InStarOnlyFunction.h"
#include <stdexcept>
#include "utils/SqrtTable.h"
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
InStarOnlyFunction::InStarOnlyFunction(string networkName, bool root) :
	NetworkAlterFunction(networkName)
{
	this->lpTable = 0;
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 * @param[simulatedState] If `true` the value() function uses the simulated
 *        state, if any or the value at the end of the period.
 */
InStarOnlyFunction::InStarOnlyFunction(string networkName, bool root, const bool simulatedState) :
	NetworkAlterFunction(networkName, simulatedState)
{
	this->lpTable = 0;
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void InStarOnlyFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpTable = this->pNetworkCache()->pInStarTable();
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void InStarOnlyFunction::initialize(const Data * pData,
	State * pState, State * pSimulatedState,  
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, pSimulatedState, period, pCache);
	this->lpTable = this->pNetworkCache()->pInStarTable();
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double InStarOnlyFunction::value(int alter) const
{
    if (this->lpTable->get(alter) > 0)
    {
        return 1.0;  // Return 1 if the value is greater than 0
    }
    return 0.0;  // Otherwise, return 0
}


/**
 * Returns the value of this function as an integer.
 */
int InStarOnlyFunction::intValue(int alter)
{
    if (this->lpTable->get(alter) > 0)
    {
        return 1;  // Return 1 if the value is greater than 0
    }
    return 0;  // Otherwise, return 0
}

}