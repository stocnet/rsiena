/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateMixedNetworkAlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * DyadicCovariateMixedNetworkAlterFunction.
 *****************************************************************************/

#include <R_ext/Print.h>
#include <stdexcept>
#include <string>

#include "DyadicCovariateMixedNetworkAlterFunction.h"
#include "data/ActorSet.h"
#include "data/DyadicCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingDyadicCovariate.h"
#include "network/Network.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "data/Data.h"
#include "MixedNetworkAlterFunction.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] firstNetworkName the name of the first network variable 
 * this function is associated with
 * @param[in] secondNetworkName the name of the second network variable 
 * this function is associated with
 * @param[in] covariateName the name of the covariate this function is
 * associated with
 */
DyadicCovariateMixedNetworkAlterFunction::DyadicCovariateMixedNetworkAlterFunction(
	string firstNetworkName, string secondNetworkName, string dyadicCovariateName) :
	MixedNetworkAlterFunction(firstNetworkName, secondNetworkName )
{
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
	this->lDyadicCovariateName = dyadicCovariateName;
}

DyadicCovariateMixedNetworkAlterFunction::~DyadicCovariateMixedNetworkAlterFunction()
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DyadicCovariateMixedNetworkAlterFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
	
	this->lpConstantCovariate =	pData->pConstantDyadicCovariate(this->lDyadicCovariateName);
	this->lpChangingCovariate =	pData->pChangingDyadicCovariate(this->lDyadicCovariateName);
	
	this->lexcludeMissings = false;
	this->lperiod = period;	

	if (this->lpConstantCovariate)
	{
		ldyadic_n = this->lpConstantCovariate->pFirstActorSet()->n();
		ldyadic_m = this->lpConstantCovariate->pSecondActorSet()->n();
	}
	else if (this->lpChangingCovariate)
	{
		ldyadic_n = this->lpChangingCovariate->pFirstActorSet()->n();
		ldyadic_m = this->lpChangingCovariate->pSecondActorSet()->n();
	}
	else
	{
		throw logic_error(
			"Dyadic covariate variable '" + this->lDyadicCovariateName + "' expected.");
	}
	
	lFirstNet_n = this->pFirstNetwork()->n();
	lFirstNet_m = this->pFirstNetwork()->m();
	lSecondNet_n = this->pSecondNetwork()->n();
	lSecondNet_m = this->pSecondNetwork()->m();	
}

/**
 * Returns the number of row nodes.
 */
int DyadicCovariateMixedNetworkAlterFunction::dyCov_n() const
{
	return ldyadic_n;
}

/**
 * Returns the number of column nodes.
 */
int DyadicCovariateMixedNetworkAlterFunction::dyCov_m() const
{
	return ldyadic_m;
}

/**
 * Returns the number of actors in the first network.
 */
int DyadicCovariateMixedNetworkAlterFunction::firstNet_n() const
{
	return lFirstNet_n;
}

/**
 * Returns the number of nodes in the second node set of the first network.
 */
int DyadicCovariateMixedNetworkAlterFunction::firstNet_m() const
{
	return lFirstNet_m;
}


/**
 * Returns the number of actors in the second network.
 */
int DyadicCovariateMixedNetworkAlterFunction::secondNet_n() const
{
	return lSecondNet_n;
}

/**
 * Returns the number of nodes in the second node set of the second network.
 */
int DyadicCovariateMixedNetworkAlterFunction::secondNet_m() const
{
	return lSecondNet_m;
}

/**
 * Returns the covariate value for the given pair of actors.
 */
double DyadicCovariateMixedNetworkAlterFunction::dyadicValue(int i, int j) const
{
	if (this->lpConstantCovariate)
	{
		return this->lpConstantCovariate->value(i, j) -
			this->lpConstantCovariate->mean();
	}
	return this->lpChangingCovariate->value(i, j, this->lperiod) -
		this->lpChangingCovariate->mean();
}

/**
 * Returns if the covariate value for the given pair of actors is missing.
 */
bool DyadicCovariateMixedNetworkAlterFunction::missing(int i, int j) const
{
	bool missing = false;

	if (this->lpConstantCovariate)
	{
		missing = this->lpConstantCovariate->missing(i, j);
	}
	else
	{
		missing = this->lpChangingCovariate->missing(i, j, this->lperiod);
	}

	return missing;
}

}
