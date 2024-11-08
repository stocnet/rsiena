/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * MixedTwoPathFunction.
 *****************************************************************************/

#include "MixedTwoStepFunction.h"
#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include "utils/SqrtTable.h"
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
MixedTwoStepFunction::MixedTwoStepFunction(std::string firstNetworkName,
                                           std::string secondNetworkName, 
                                           Direction firstDirection, 
                                           Direction secondDirection, double par) :
	MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
{
	this->lpTable = 0;
	this->ldirection1 = firstDirection;
	this->ldirection2 = secondDirection;
	this->lroot = (fabs(par-2) < 0.001);
	this->ltrunc = (fabs(par-3) < 0.001);
	this->lsqrtTable = SqrtTable::instance();
}

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
MixedTwoStepFunction::MixedTwoStepFunction(std::string firstNetworkName,
                                           std::string secondNetworkName, 
                                           Direction firstDirection, 
                                           Direction secondDirection, 
										   double par, bool simulatedState) :
	MixedNetworkAlterFunction(firstNetworkName, secondNetworkName, simulatedState)
{
	this->lpTable = 0;
	this->ldirection1 = firstDirection;
	this->ldirection2 = secondDirection;
	this->lroot = (fabs(par-2) < 0.001);
	this->ltrunc = (fabs(par-3) < 0.001);
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void MixedTwoStepFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
	if (ldirection1 == EITHER && ldirection2 == EITHER)
		this->lpTable = this->pTwoNetworkCache()->pEETable();
	if (ldirection1 == FORWARD && ldirection2 == EITHER)
		this->lpTable = this->pTwoNetworkCache()->pFETable();
	if (ldirection1 == FORWARD && ldirection2 == RECIPROCAL)
		this->lpTable = this->pTwoNetworkCache()->pFRTable();
	if (ldirection1 == EITHER && ldirection2 == RECIPROCAL)
		this->lpTable = this->pTwoNetworkCache()->pERTable();
	if (ldirection1 == FORWARD && ldirection2 == FORWARD)
		this->lpTable = this->pTwoNetworkCache()->pTwoPathTable();
	if (ldirection1 == BACKWARD && ldirection2 == FORWARD)
		this->lpTable = this->pTwoNetworkCache()->pOutStarTable();
	if (ldirection1 == FORWARD && ldirection2 == BACKWARD)
		this->lpTable = this->pTwoNetworkCache()->pInStarTable();
	if (ldirection1 == RECIPROCAL && ldirection2 == FORWARD)
		this->lpTable = this->pTwoNetworkCache()->pRFTable();

	if(this->lpTable == 0)
		throw std::invalid_argument( "MixedTwoStepFunction expects different direction parameters" );
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] pSimulatedState the current simulated state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void MixedTwoStepFunction::initialize(const Data * pData,
	State * pState, State * pSimulatedState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, pSimulatedState, period, pCache);
	if (ldirection1 == EITHER && ldirection2 == EITHER)
		this->lpTable = this->pTwoNetworkCache()->pEETable();
	if (ldirection1 == FORWARD && ldirection2 == EITHER)
		this->lpTable = this->pTwoNetworkCache()->pFETable();
	if (ldirection1 == FORWARD && ldirection2 == RECIPROCAL)
		this->lpTable = this->pTwoNetworkCache()->pFRTable();
	if (ldirection1 == EITHER && ldirection2 == RECIPROCAL)
		this->lpTable = this->pTwoNetworkCache()->pERTable();
	if (ldirection1 == FORWARD && ldirection2 == FORWARD)
		this->lpTable = this->pTwoNetworkCache()->pTwoPathTable();
	if (ldirection1 == BACKWARD && ldirection2 == FORWARD)
		this->lpTable = this->pTwoNetworkCache()->pOutStarTable();
	if (ldirection1 == FORWARD && ldirection2 == BACKWARD)
		this->lpTable = this->pTwoNetworkCache()->pInStarTable();
	if (ldirection1 == RECIPROCAL && ldirection2 == FORWARD)
		this->lpTable = this->pTwoNetworkCache()->pRFTable();

	if(this->lpTable == 0)
		throw std::invalid_argument( "MixedTwoStepFunction expects different direction parameters" );
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double MixedTwoStepFunction::value(int alter) const
{
	if (this->lroot)
	{
		return this->lsqrtTable->sqrt(this->lpTable->get(alter));
	}
	else if (this->ltrunc)
	{
		if (this->lpTable->get(alter) >= 1)
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


/**
 * Returns the value of this function as an integer.
 */
int MixedTwoStepFunction::intValue(int alter)
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
