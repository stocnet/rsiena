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

	if (!this->lpConstantCovariate && !this->lpChangingCovariate)
	{
		throw logic_error(
			"Dyadic covariate variable '" + this->lDyadicCovariateName + "' expected.");
	}
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
