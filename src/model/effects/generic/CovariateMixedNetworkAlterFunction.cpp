/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateMixedNetworkAlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateMixedNetworkAlterFunction.
 *****************************************************************************/

#include <R_ext/Print.h>
#include <stdexcept>
#include <string>

#include "CovariateMixedNetworkAlterFunction.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/variables/BehaviorVariable.h"
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
CovariateMixedNetworkAlterFunction::CovariateMixedNetworkAlterFunction(
	string firstNetworkName, string secondNetworkName, string covariateName) :
	MixedNetworkAlterFunction(firstNetworkName, secondNetworkName )
{
	this->lcovariateName = covariateName;
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
	this->lpBehaviorData = 0;
	this->lvalues = 0;
}

CovariateMixedNetworkAlterFunction::~CovariateMixedNetworkAlterFunction()
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateMixedNetworkAlterFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpConstantCovariate = pData->pConstantCovariate(this->lcovariateName);
	this->lpChangingCovariate = pData->pChangingCovariate(this->lcovariateName);
	this->lpBehaviorData = pData->pBehaviorData(this->lcovariateName);
	this->lvalues = pState->behaviorValues(this->lcovariateName);
	this->lperiod = period;

	if (!this->lpConstantCovariate && !this->lpChangingCovariate
			&& !(this->lpBehaviorData && this->lvalues))
	{
		throw logic_error("Covariate or dependent behavior variable '" +
				this->lcovariateName + "' expected.");
	}
}

/**
 * Returns the covariate value for the given actor.
 */
double CovariateMixedNetworkAlterFunction::value(int alter) const
{
	if (this->lpConstantCovariate)
	{
		return this->lpConstantCovariate->value(alter);
	}
	if (this->lpChangingCovariate)
	{
		return this->lpChangingCovariate->value(alter, this->lperiod);
	}
	return this->lvalues[alter] - this->lpBehaviorData->overallMean();
}

/**
 * Returns if the covariate value for the given actor is missing.
 */
bool CovariateMixedNetworkAlterFunction::missing(int i) const
{
	bool missing = false;

	if (this->lpConstantCovariate)
	{
		missing = this->lpConstantCovariate->missing(i);
	}
	else if (this->lpChangingCovariate)
	{
		missing = this->lpChangingCovariate->missing(i, this->lperiod);
	}
	else
	{
		missing = this->lpBehaviorData->missing(this->lperiod, i);
	}

	return missing;
}

/**
 * Returns the centered similarity of the given actors.
 */
double CovariateMixedNetworkAlterFunction::similarity(int i, int j) const
{
	double similarity = 0;

	if (this->lpConstantCovariate)
	{
		similarity =
			this->lpConstantCovariate->similarity(
				this->lpConstantCovariate->value(i),
				this->lpConstantCovariate->value(j));
	}
	else if (this->lpChangingCovariate)
	{
		similarity =
			this->lpChangingCovariate->similarity(this->value(i),
				this->value(j));
	}
	else
	{
		similarity =
			this->lpBehaviorData->similarity(this->lvalues[i],
				this->lvalues[j]);
	}
	return similarity;
}

/**
 * Returns the constant covariate associated with this effect.
 */
ConstantCovariate * CovariateMixedNetworkAlterFunction::pConstantCovariate() const
{
		return this->lpConstantCovariate;
}

/**
 * Returns the changing covariate associated with this effect.
 */
ChangingCovariate * CovariateMixedNetworkAlterFunction::pChangingCovariate() const
{
		return this->lpChangingCovariate;
}

/**
 * Returns the changing covariate associated with this effect.
 */
BehaviorLongitudinalData * CovariateMixedNetworkAlterFunction::pBehaviorData() const
{
		return this->lpBehaviorData;
}

}
