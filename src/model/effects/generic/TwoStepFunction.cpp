/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TwoStepFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * TwoStepFunction.
 *****************************************************************************/

#include "TwoStepFunction.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include <stdexcept>

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
TwoStepFunction::TwoStepFunction(std::string networkName, Direction direction1, Direction direction2) :
	OneModeNetworkAlterFunction(networkName)
{
	this->lpTable = 0;
	this->ldirection1 = direction1;
	this->ldirection2 = direction2;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void TwoStepFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	OneModeNetworkAlterFunction::initialize(pData, pState, period, pCache);

	if (ldirection1 == FORWARD && ldirection2 == RECIPROCAL)
		this->lpTable = this->pNetworkCache()->pFRTable();
	if (ldirection1 == FORWARD && ldirection2 == FORWARD)
		this->lpTable = this->pNetworkCache()->pTwoPathTable();


	if(this->lpTable == 0)
		throw std::invalid_argument( "TwoStepFunction expects different direction parameters" );
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double TwoStepFunction::value(int alter) const
{
	return this->lpTable->get(alter);
}


/**
 * Returns the value of this function as an integer.
 */
int TwoStepFunction::intValue(int alter)
{
	return this->lpTable->get(alter);
}

}
