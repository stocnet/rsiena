/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateNetworkFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateNetworkAlterFunction.
 *****************************************************************************/

#include <stdexcept>
#include <string>
#include <cmath>

#include "CovariateNetworkAlterFunction.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/variables/BehaviorVariable.h"
#include "network/Network.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "data/Data.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 * @param[in] covariateName the name of the covariate this function is
 * associated with
 */
CovariateNetworkAlterFunction::CovariateNetworkAlterFunction(string networkName,
string covariateName) :
	NetworkAlterFunction(networkName)
{
	this->lcovariateName = covariateName;
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
	this->lpBehaviorData = 0;
	this->lvalues = 0;
}

CovariateNetworkAlterFunction::~CovariateNetworkAlterFunction()
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateNetworkAlterFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);

	this->lpConstantCovariate = pData->pConstantCovariate(this->lcovariateName);
	this->lpChangingCovariate = pData->pChangingCovariate(this->lcovariateName);
	this->lpBehaviorData = pData->pBehaviorData(this->lcovariateName);
	this->lvalues = pState->behaviorValues(this->lcovariateName);

	this->lperiod = period;

	if (!this->lpConstantCovariate &&
		!this->lpChangingCovariate &&
		!(this->lpBehaviorData && this->lvalues))
	{
		throw logic_error("Covariate or dependent behavior variable '" +
			this->lcovariateName +
			"' expected.");
	}
}

/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling CovariateNetworkAlterFunction::value(...).
 */
void CovariateNetworkAlterFunction::preprocessEgo(int ego)
{
	NetworkAlterFunction::preprocessEgo(ego);
}

/**
 * Returns the overall mean of covvalue
 */
double CovariateNetworkAlterFunction::covmean() const
{
	double themean = 0;

	if (this->lpConstantCovariate)
	{
		themean = this->lpConstantCovariate->mean();
	}
	else if (this->lpChangingCovariate)
	{
		themean = this->lpChangingCovariate->mean();
	}
	// else lpBehaviorData: values are already centered,
	 // see CovariateNetworkAlterFunction::covvalue
	 // themean is already 0

	return themean;
}

/**
 * Returns the covariate value for the given actor.
 * For behavior, returns the centered value.
 */
double CovariateNetworkAlterFunction::covvalue(int alter) const
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
 * Returns the covariate value for the given actor, rounded to integer.
 * For behavior, this is the non-centered value.
 */
int CovariateNetworkAlterFunction::covIntValue(int i) const
{
	if (this->lpBehaviorData)
	{
		return this->lvalues[i] ;
	}
	else
	{
		return int(round(this->covvalue(i)));
	}
}


/**
 * Returns if the covariate value for the given actor is missing.
 */
bool CovariateNetworkAlterFunction::missing(int i) const
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
 * Returns the covariate number of cases.
 */
int CovariateNetworkAlterFunction::covariateN() const
{
	int ncov = 0;

	if (this->lpConstantCovariate)
	{
		ncov = this->lpConstantCovariate->covariateN();

	}
	else if (this->lpChangingCovariate)
	{
		ncov = this->lpChangingCovariate->covariateN();
	}
	else
	{
		ncov = this->lpBehaviorData->observationCount();
	}
	
	return ncov;
}

/**
 * Returns the first covariate minimum value.
 * For behavior, this is the minimum non-centered value.
 */
double CovariateNetworkAlterFunction::covariateMinimum() const
{
	double mini = 0;

	if (this->lpConstantCovariate)
	{
		mini = this->lpConstantCovariate->min();

	}
	else if (this->lpChangingCovariate)
	{
		mini = this->lpChangingCovariate->min();
	}
	else
	{
		mini = this->lpBehaviorData->min();
	}
	
	return mini;
}


/**
 * Returns the first covariate maximum value.
 * For behavior, this is the maximum non-centered value.
 */
double CovariateNetworkAlterFunction::covariateMaximum() const
{
	double maxi = 0;

	if (this->lpConstantCovariate)
	{
		maxi = this->lpConstantCovariate->max();

	}
	else if (this->lpChangingCovariate)
	{
		maxi = this->lpChangingCovariate->max();
	}
	else
	{
		maxi = this->lpBehaviorData->max();
	}
	
	return maxi;
}


/**
 * Returns the centered similarity of the given actors.
 */
double CovariateNetworkAlterFunction::actor_similarity(int i, int j) const
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
			this->lpChangingCovariate->similarity(this->lpChangingCovariate->value(i, this->lperiod),
				this->lpChangingCovariate->value(j, this->lperiod));
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
ConstantCovariate * CovariateNetworkAlterFunction::pConstantCovariate() const
{
	return this->lpConstantCovariate;
}

/**
 * Returns the changing covariate associated with this effect.
 */
ChangingCovariate * CovariateNetworkAlterFunction::pChangingCovariate() const
{
	return this->lpChangingCovariate;
}

/**
 * Returns the changing covariate associated with this effect.
 */
BehaviorLongitudinalData * CovariateNetworkAlterFunction::pBehaviorData() const
{
	return this->lpBehaviorData;
}

}
