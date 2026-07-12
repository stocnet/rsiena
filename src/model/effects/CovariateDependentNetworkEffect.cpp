/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDependentNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDependentNetworkEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "CovariateDependentNetworkEffect.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ContinuousLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/ContinuousVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
CovariateDependentNetworkEffect::CovariateDependentNetworkEffect(
		const EffectInfo * pEffectInfo) :
		NetworkEffect(pEffectInfo), //
		lSimulatedOffset(0), //
		lpConstantCovariate(0), //
		lpChangingCovariate(0), //
		lpBehaviorData(0), //
        lpContinuousData(0),
		lvalues(0),
        lcontinuousValues(0) {
}

/**
 * Constructor.
 *
 * @param pEffectInfo The effect info.
 * in addition, for gmom:
 * @param simulatedState If `true` the value(), missing() and actor_similarity()
 *        functions uses the simulated state, if any or the value at the end
 *        of the period.
 */
CovariateDependentNetworkEffect::CovariateDependentNetworkEffect(
		const EffectInfo * pEffectInfo, bool simulatedState) :
		NetworkEffect(pEffectInfo), //
		lSimulatedOffset(simulatedState ? 1 : 0), //
		lpConstantCovariate(0), //
		lpChangingCovariate(0), //
		lpBehaviorData(0), //
        lpContinuousData(0),
		lvalues(0),
        lcontinuousValues(0) {
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the state of the dependent variables at the beginning of the period
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateDependentNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
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
			name + "' expected.");
	}
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the state of the dependent variables at the beginning of the period
 * @param[in] pSimulatedState the current simulated state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateDependentNetworkEffect::initialize(const Data* pData,
		State* pState, State* pSimulatedState, int period, Cache* pCache) {
	NetworkEffect::initialize(pData, pState, pSimulatedState, period, pCache);
	const string name = this->pEffectInfo()->interactionName1();

	// For covariates lSimulatedOffset is applied in values().
	this->lpConstantCovariate = pData->pConstantCovariate(name);
	this->lpChangingCovariate = pData->pChangingCovariate(name);

	this->lpBehaviorData = pData->pBehaviorData(name);
    this->lpContinuousData = pData->pContinuousData(name);
	// If the covariate is a behaviour take the right state.
	if (lSimulatedOffset == 1) {
		this->lvalues = pSimulatedState->behaviorValues(name);
        this->lcontinuousValues = pSimulatedState->continuousValues(name);
	} else {
		this->lvalues = pState->behaviorValues(name);
        this->lcontinuousValues = pState->continuousValues(name);
	}

	if (!this->lpConstantCovariate &&
        !this->lpChangingCovariate &&
        !(this->lpBehaviorData && this->lvalues) &&
        !(this->lpContinuousData && this->lcontinuousValues)) {
		throw logic_error(
				"Covariate or dependent behavior variable '" + name + "' expected.");
	}
}

/**
 * Returns the covariate value for the given actor.
 * Note that this differs from what is done in BehaviorEffect,
 * where value is the non-centered and centeredValue is the centered value.
 */
double CovariateDependentNetworkEffect::value(const int i) const
{
    double value = 0;
    
    if (this->lpConstantCovariate)
    {
		value = this->lpConstantCovariate->value(i);
	}
	else if (this->lpChangingCovariate)
    {
		value = this->lpChangingCovariate->value(i, this->period() + lSimulatedOffset);
	}
    else if (this->lpBehaviorData)
    {
        value = this->lvalues[i] - this->lpBehaviorData->overallMean();
    }
    else
    {
        value = this->lcontinuousValues[i] - this->lpContinuousData->overallMean();
    }
    
    return value;
}


/**
 * Returns the non-centered (raw) covariate value for the given actor.
 * For a constant or changing covariate that was declared centered, the stored
 * value has the observed mean subtracted, so it is added back; a non-centered
 * covariate is already raw and is returned unchanged.
 * For behavior variables this returns lvalues[alter] without mean subtraction.
*/
double CovariateDependentNetworkEffect::uncenteredValue(const int i) const
{
	double value = 0;

	if (this->lpConstantCovariate)
	{
		value = this->lpConstantCovariate->value(i);
		if (this->lpConstantCovariate->centered())
		{
			value += this->lpConstantCovariate->obsmean();
		}
	}
	else if (this->lpChangingCovariate)
	{
		value = this->lpChangingCovariate->value(i, this->period() + lSimulatedOffset);
		if (this->lpChangingCovariate->centered())
		{
			value += this->lpChangingCovariate->obsmean();
		}
	}
	else if (this->lpBehaviorData)
	{
		value = this->lvalues[i];
	}
	else
	{
		value = this->lcontinuousValues[i];
	}

	return value;
}

/**
 * Returns the covariate/behavior value of a given actor in the requested
 * form, e.g. centered around the overall mean of all observed values or
 * not. Mirrors BehaviorEffect::resolvedValue.
 * @param[in] i the index of the actor
 * @param[in] nc whether the value should not be centered
 */
double CovariateDependentNetworkEffect::resolvedValue(const int i, bool nc) const
{
	return nc ? this->uncenteredValue(i) : this->value(i);
}

/**
 * Returns if the covariate value for the given actor is missing.
 */
bool CovariateDependentNetworkEffect::missing(int i) const
{
    bool missing = false;

    if (this->lpConstantCovariate)
    {
		missing = this->lpConstantCovariate->missing(i);
	}
	else if (this->lpChangingCovariate)
    {
		missing = this->lpChangingCovariate->missing(i, this->period() + lSimulatedOffset);
	}
    else if (this->lpBehaviorData)
    {
        missing = this->lpBehaviorData->missing(this->period() + lSimulatedOffset, i);
        // This means: if simulated state, period+1; else period.
    }
    else
    {
        missing = this->lpContinuousData->missing(this->period() + lSimulatedOffset, i);

    }
    
    return missing;
}


/**
 * Returns the covariate number of cases.
 * For behavior, this is the minimum non-centered value.
 */
int CovariateDependentNetworkEffect::covarN() const
{
	int covN = 0;
	if (this->lpConstantCovariate)
	{
		covN = this->lpConstantCovariate->covariateN();
	}
	else if (this->lpChangingCovariate)
	{
		covN = this->lpChangingCovariate->covariateN();
	}
	else
	{
		covN = this->lpBehaviorData->observationCount();
	}
	return covN;
}

/**
 * Returns the covariate minimum value.
 * For behavior, this is the minimum non-centered value.
 */
double CovariateDependentNetworkEffect::covariateMinimum() const
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
 * Returns the covariate maximum value.
 * For behavior, this is the maximum non-centered value.
 */
double CovariateDependentNetworkEffect::covariateMaximum() const
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
double CovariateDependentNetworkEffect::actor_similarity(int i, int j) const
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
	else if (this->lpBehaviorData)
	{
		similarity =
			this->lpBehaviorData->similarity(this->lvalues[i],
				this->lvalues[j]);
	}
    else
    {
        similarity =
            this->lpContinuousData->similarity(this->lcontinuousValues[i],
                this->lcontinuousValues[j]);
    }

	return similarity;
}

/**
 * Returns the constant covariate associated with this effect.
 */
ConstantCovariate * CovariateDependentNetworkEffect::pConstantCovariate() const
{
		return this->lpConstantCovariate;

}

/**
 * Returns the changing covariate associated with this effect.
 */
ChangingCovariate * CovariateDependentNetworkEffect::pChangingCovariate() const
{
		return this->lpChangingCovariate;

}

/**
 * Returns the behavioral variable associated with this effect.
 */
BehaviorLongitudinalData * CovariateDependentNetworkEffect::pBehaviorData() const
{
		return this->lpBehaviorData;
}

/**
 * Returns the continuous behavioral variable associated with this effect.
 */
ContinuousLongitudinalData * CovariateDependentNetworkEffect::pContinuousData() const
{
        return this->lpContinuousData;
}

}
