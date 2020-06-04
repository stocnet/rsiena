/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDependentBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDependentBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "CovariateDependentBehaviorEffect.h"
#include "data/Data.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
CovariateDependentBehaviorEffect::CovariateDependentBehaviorEffect(
	const EffectInfo * pEffectInfo) : BehaviorEffect(pEffectInfo)
{
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
	this->lpBehaviorData = 0;
	this->linteractionValues = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateDependentBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	BehaviorEffect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->interactionName1();

	this->lpConstantCovariate = pData->pConstantCovariate(name);
	this->lpChangingCovariate = pData->pChangingCovariate(name);
	this->lpBehaviorData = pData->pBehaviorData(name);
	this->linteractionValues = pState->behaviorValues(name);

	if (!this->lpConstantCovariate &&
		!this->lpChangingCovariate &&
		!(this->lpBehaviorData && this->linteractionValues))
	{
		throw logic_error("Covariate or dependent behavior variable '" +
			name +
			"' expected.");
	}
}


/**
 * Returns the observed overall mean of covariateValue
 */
double CovariateDependentBehaviorEffect::covariateMean() const
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
	 // see CovariateNetworkAlterFunction::value
	 // themean is already 0

	return themean;
}

/**
 * Returns the covariate value for the given actor.
 */
double CovariateDependentBehaviorEffect::covariateValue(int i) const
{
	double value = 0;

	if (this->lpConstantCovariate)
	{
		value = this->lpConstantCovariate->value(i);
	}
	else if (this->lpChangingCovariate)
	{
		value = this->lpChangingCovariate->value(i, this->period());
	}
	else
	{
		value = this->linteractionValues[i] -
			this->lpBehaviorData->overallMean();
	}

	return value;
}


/**
 * Returns if the covariate value for the given actor is missing at the
 * given observation.
 */
bool CovariateDependentBehaviorEffect::missingCovariate(int i,
	int observation) const
{
	bool missing = false;

	if (this->lpConstantCovariate)
	{
		missing = this->lpConstantCovariate->missing(i);
	}
	else if (this->lpChangingCovariate)
	{
		missing = this->lpChangingCovariate->missing(i, observation);
	}
	else
	{
		missing = this->lpBehaviorData->missing(observation, i);
	}

	return missing;
}

/**
 * Returns if the covariate value for the given actor is missing at the
 * given observation or the next one. The next one is only used for
 * behavior variable. It may not exist for changing covariate so is never used.
 */
bool CovariateDependentBehaviorEffect::missingCovariateEitherEnd(int i,
	int observation) const
{
	bool missing = false;

	if (this->lpConstantCovariate)
	{
		missing = this->lpConstantCovariate->missing(i);
	}
	else if (this->lpChangingCovariate)
	{
		missing = this->lpChangingCovariate->missing(i, observation);
	}
	else
	{
		missing = this->lpBehaviorData->missing(observation, i) ||
		this->lpBehaviorData->missing(observation + 1, i);
	}

	return missing;
}
}
