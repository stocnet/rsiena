/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateDependentNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DyadicCovariateDependentNetworkEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <R_ext/Print.h>
#include "DyadicCovariateDependentNetworkEffect.h"
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
DyadicCovariateDependentNetworkEffect::DyadicCovariateDependentNetworkEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DyadicCovariateDependentNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->interactionName1();

	this->lpConstantCovariate =	pData->pConstantDyadicCovariate(name);
	this->lpChangingCovariate =	pData->pChangingDyadicCovariate(name);

	this->lexcludeMissings = false;

	if (!this->lpConstantCovariate && !this->lpChangingCovariate)
	{
		throw logic_error(
			"Dyadic covariate variable '" + name + "' expected.");
	}
}


/**
 * Returns the covariate value for the given pair of actors.
 */
double DyadicCovariateDependentNetworkEffect::value(int i, int j) const
{
	if (this->lpConstantCovariate)
	{
		return this->lpConstantCovariate->value(i, j) -
			this->lpConstantCovariate->mean();
	}
	return this->lpChangingCovariate->value(i, j, this->period()) -
		this->lpChangingCovariate->mean();
}


/**
 * Returns if the covariate value for the given pair of actors is missing.
 */
bool DyadicCovariateDependentNetworkEffect::missing(int i, int j) const
{
	bool missing = false;

	if (this->lpConstantCovariate)
	{
		missing = this->lpConstantCovariate->missing(i, j);
	}
	else
	{
		missing = this->lpChangingCovariate->missing(i, j, this->period());
	}

	return missing;
}

/**
 * Returns if the associated covariate is a constant covariate or not
 */
bool DyadicCovariateDependentNetworkEffect::constantCovariate() const
{

	if (this->lpConstantCovariate)
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
	DyadicCovariateDependentNetworkEffect::rowValues(int i) const
{
	if (this->lpConstantCovariate)
	{
		return this->lpConstantCovariate->rowValues(i);
	}
	else
	{
		//	Rprintf("%d %d effect \n", i, this->lexcludeMissings);
		return this->lpChangingCovariate->rowValues(i, this->period(),
			this->lexcludeMissings);
	}
}


/**
 * Returns an iterator over non-zero non-missing values of the given column
 * of the covariate.
 */
DyadicCovariateValueIterator
	DyadicCovariateDependentNetworkEffect::columnValues(int j) const
{
	if (this->lpConstantCovariate)
	{
		return this->lpConstantCovariate->columnValues(j);
	}
	else
	{
		//	Rprintf("%d effect \n", this->lexcludeMissings);
		return this->lpChangingCovariate->columnValues(j, this->period(),
			this->lexcludeMissings);
	}
}

/**
 * This method is called at the start of the calculation of the statistic.
 */
void DyadicCovariateDependentNetworkEffect::initializeStatisticCalculation()
{
		this->lexcludeMissings = true;
}
/**
 * This method is called at the end of the calculation of the statistic.
 */
void DyadicCovariateDependentNetworkEffect::cleanupStatisticCalculation()
{
		this->lexcludeMissings = false;
}


}
