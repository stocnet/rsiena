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
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "model/variables/BehaviorVariable.h"


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
		lvalues(0) {
}

/**
 * Constructor.
 *
 * @param pEffectInfo The effect info.
 * @param simulatedState If `true` the value(), missing() and similarity()
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
		lvalues(0) {
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
	this->lvalues = pState->behaviorValues(name);

	if (!this->lpConstantCovariate &&
		!this->lpChangingCovariate &&
		!(this->lpBehaviorData && this->lvalues))
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
	// If the covariate is a behaviour take the right state.
	if (lSimulatedOffset == 1) {
		this->lvalues = pSimulatedState->behaviorValues(name);
	} else {
		this->lvalues = pState->behaviorValues(name);
	}

	if (!this->lpConstantCovariate && !this->lpChangingCovariate
			&& !(this->lpBehaviorData && this->lvalues)) {
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
	if (this->lpConstantCovariate) {
		return this->lpConstantCovariate->value(i);
	}
	if (this->lpChangingCovariate) {
		return this->lpChangingCovariate->value(i, this->period() + lSimulatedOffset);
	}
  return this->lvalues[i] - this->lpBehaviorData->overallMean();
}


/**
 * Returns if the covariate value for the given actor is missing.
 */
bool CovariateDependentNetworkEffect::missing(int i) const
{
	if (this->lpConstantCovariate) {
		return this->lpConstantCovariate->missing(i);
	}
	if (this->lpChangingCovariate) {
		return this->lpChangingCovariate->missing(i, this->period() + lSimulatedOffset);
	}
  return this->lpBehaviorData->missing(this->period() + lSimulatedOffset, i);
// This means: if simulated state, period+1; else period.
}


/**
 * Returns the centered similarity of the given actors.
 */
double CovariateDependentNetworkEffect::similarity(int i, int j) const
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

}
