/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InStarFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * InStarFunction.
 *****************************************************************************/

#include "InStarFunction.h"
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
InStarFunction::InStarFunction(string networkName, bool root, bool only) :
	NetworkAlterFunction(networkName)
{
	this->lpTable = 0;
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
	this->lonly = only;
}

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 * @param[simulatedState] If `true` the value() function uses the simulated
 *        state, if any or the value at the end of the period.
 */
InStarFunction::InStarFunction(string networkName, bool root, bool only, 
							const bool simulatedState) :
	NetworkAlterFunction(networkName, simulatedState)
{
	this->lpTable = 0;
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
	this->lonly = only;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void InStarFunction::initialize(const Data * pData,
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
void InStarFunction::initialize(const Data * pData,
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
double InStarFunction::value(int alter) const
{
	if (this->lroot)
	{
		return this->lsqrtTable->sqrt(this->lpTable->get(alter));
	}
	else if (this->lonly)
	{
		if (this->lpTable->get(alter) > 0)
		{
			return 1.0;
		}
		else
		{
			return 0.0;
		}
	}
	else
	{
		return this->lpTable->get(alter);
	}
}


/**
 * Returns the value of this function as an integer.
 */
int InStarFunction::intValue(int alter)
{
	if (this->lroot)
	{
		throw logic_error("Square roots are not integer values");
	}
	else if (this->lonly)
	{
		if (this->lpTable->get(alter) > 0)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return this->lpTable->get(alter);
	}
}

}
