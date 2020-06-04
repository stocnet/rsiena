/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Behavior Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateAndNetworkBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DyadicCovariateAndNetworkBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <R_ext/Print.h>
#include "DyadicCovariateAndNetworkBehaviorEffect.h"
#include "data/Data.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/variables/BehaviorVariable.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/DyadicCovariateValueIterator.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
DyadicCovariateAndNetworkBehaviorEffect::DyadicCovariateAndNetworkBehaviorEffect(
	const EffectInfo * pEffectInfo) : NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->lpConstantDyadicCovariate = 0;
	this->lpChangingDyadicCovariate = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DyadicCovariateAndNetworkBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkDependentBehaviorEffect::initialize(pData, pState, period, pCache);
	string name1 = this->pEffectInfo()->interactionName1();
	string name2 = this->pEffectInfo()->interactionName2();

	this->lpConstantDyadicCovariate = pData->pConstantDyadicCovariate(name2);
	this->lpChangingDyadicCovariate = pData->pChangingDyadicCovariate(name2);
	this->lpBehaviorData = pData->pBehaviorData(name1);
	this->lexcludeMissings = false;

	if (!this->lpConstantDyadicCovariate && !this->lpChangingDyadicCovariate)
	{
		throw logic_error(
			"Dyadic covariate variable '" + name2 + "' expected.");
	}
}

/**
 * Returns the covariate value for the given pair of actors.
 */
double DyadicCovariateAndNetworkBehaviorEffect::dycoValue(int i, int j) const
{
	double value = 0;

	if (this->lpConstantDyadicCovariate)
	{
		value = this->lpConstantDyadicCovariate->value(i, j) -
			this->lpConstantDyadicCovariate->mean();
	}
	else
	{
		value = this->lpChangingDyadicCovariate->value(i, j, this->period()) -
			this->lpChangingDyadicCovariate->mean();
	}

	return value;
}


/**
 * Returns if the covariate value for the given pair of actors is missing.
 */
bool DyadicCovariateAndNetworkBehaviorEffect::missingDyCo(int i, int j) const
{
	bool missing = false;

	if (this->lpConstantDyadicCovariate)
	{
		missing = this->lpConstantDyadicCovariate->missing(i, j);
	}
	else
	{
		missing = this->lpChangingDyadicCovariate->missing(i, j, this->period());
	}

	return missing;
}


/**
 * Returns if the associated covariate is a constant covariate or not
 */
bool DyadicCovariateAndNetworkBehaviorEffect::constantDyadicCovariate() const
{
	if (this->lpConstantDyadicCovariate)
	{
		return true;
	}
	else
	{
		return false;
	}
}


/**
 * Returns an iterator over non-zero non-missing values of the given row
 * of the covariate.
 */
DyadicCovariateValueIterator
	DyadicCovariateAndNetworkBehaviorEffect::rowValues(int i) const
{
	if (this->lpConstantDyadicCovariate)
	{
		return this->lpConstantDyadicCovariate->rowValues(i);
	}
	else
	{
		//	Rprintf("%d %d effect \n", i, this->lexcludeMissings);
		return this->lpChangingDyadicCovariate->rowValues(i, this->period(),
			this->lexcludeMissings);
	}
}

/**
 * Returns an iterator over non-zero non-missing values of the given column
 * of the covariate.
 */
DyadicCovariateValueIterator
	DyadicCovariateAndNetworkBehaviorEffect::columnValues(int j) const
{
	if (this->lpConstantDyadicCovariate)
	{
		return this->lpConstantDyadicCovariate->columnValues(j);
	}
	else
	{
		//	Rprintf("%d effect \n", this->lexcludeMissings);
		return this->lpChangingDyadicCovariate->columnValues(j, this->period(),
			this->lexcludeMissings);
	}
}

/**
 * This method is called at the start of the calculation of the
 * evaluationStatistic, endowmentStatistic, and creationStatistic
 */
void DyadicCovariateAndNetworkBehaviorEffect::initializeStatisticCalculation()
{
		this->lexcludeMissings = true;
// Prevents having to check missingness in egoStatistic()
}

/**
 * This method is called at the end of the calculation of the
 * evaluationStatistic, endowmentStatistic, and creationStatistic
 */
void DyadicCovariateAndNetworkBehaviorEffect::cleanupStatisticCalculation()
{
		this->lexcludeMissings = false;
}
}
