/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDependentContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDependentContinuousEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "CovariateDependentContinuousEffect.h"
#include "data/Data.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ContinuousLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/ContinuousVariable.h"

namespace siena
{

/**
 * Constructor.
 */
CovariateDependentContinuousEffect::CovariateDependentContinuousEffect(
	const EffectInfo * pEffectInfo) : ContinuousEffect(pEffectInfo)
{
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
	this->lpBehaviorData = 0;
	this->lpContinuousData = 0;
	this->lvalues = 0;
	this->lcontinuousValues = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateDependentContinuousEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	ContinuousEffect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->interactionName1();

	this->lpConstantCovariate = pData->pConstantCovariate(name);
	this->lpChangingCovariate = pData->pChangingCovariate(name);
	this->lpBehaviorData = pData->pBehaviorData(name);
	this->lpContinuousData = pData->pContinuousData(name);
	this->lvalues = pState->behaviorValues(name);
	this->lcontinuousValues = pState->continuousValues(name);

	
	if (!this->lpConstantCovariate &&
		!this->lpChangingCovariate &&
		!(this->lpBehaviorData && this->lvalues) &&
		!(this->lpContinuousData && this->lcontinuousValues))
	{
		throw logic_error("Covariate or dependent behavior variable '" +
			name +
			"' expected.");
	}
}


/**
 * Returns the covariate value for the given actor.
 */
double CovariateDependentContinuousEffect::covariateValue(int i) const
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
	else if (this->lpBehaviorData)
	{
		value = this->lvalues[i] - this->lpBehaviorData->overallMean();
	} 
	else // no extra centering: the interaction value should be 'as centered'
		 // as the data value itself (conform SDE)
	{
		value = this->lcontinuousValues[i];
	}

	return value;
}


/**
 * Returns if the covariate value for the given actor is missing at the
 * given observation.
 */
bool CovariateDependentContinuousEffect::missingCovariate(int i,
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
	else if (this->lpBehaviorData)
	{
		missing = this->lpBehaviorData->missing(observation, i);
	}
	else
	{
		missing = this->lpContinuousData->missing(observation, i);
	}
	
	return missing;
}

/**
 * Returns if the covariate value for the given actor is missing at the
 * given observation or the next one. The next one is only used for
 * behavior variables (continuous or discrete). It may not exist for 
 * changing covariate so is never used.
 */
bool CovariateDependentContinuousEffect::missingCovariateEitherEnd(int i,
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
	else if (this->lpBehaviorData)
	{
		missing = this->lpBehaviorData->missing(observation, i) ||
		this->lpBehaviorData->missing(observation + 1, i);
	}
	else
	{
		missing = this->lpContinuousData->missing(observation, i) ||
		this->lpContinuousData->missing(observation + 1, i);
	}
	return missing;
}
}
