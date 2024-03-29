/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocatedTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * ReciprocatedTwoPathFunction.
 *****************************************************************************/

#include "ReciprocatedTwoPathFunction.h"
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
ReciprocatedTwoPathFunction::ReciprocatedTwoPathFunction(string networkName,
															bool root) :
	OneModeNetworkAlterFunction(networkName)
{
	this->lpTable = 0;;
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
void ReciprocatedTwoPathFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	OneModeNetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpTable = this->pNetworkCache()->pRRTable();
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double ReciprocatedTwoPathFunction::value(int alter) const
{
	if (this->lroot)
	{
		return this->lsqrtTable->sqrt(this->lpTable->get(alter));
	}
	else
	{
		return this->lpTable->get(alter);
	}
}


/**
 * Returns the value of this function as an integer.
 */
int ReciprocatedTwoPathFunction::intValue(int alter)
{
	if (this->lroot)
	{
		throw logic_error("Square roots are not integer values");
	}
	else
	{
		return this->lpTable->get(alter);
	}
}

}
