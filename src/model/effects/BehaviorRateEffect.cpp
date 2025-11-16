/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorRateEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BehaviorRateEffect.
 *****************************************************************************/

#include <stdexcept>

#include "BehaviorRateEffect.h"
#include "data/Data.h"
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
BehaviorRateEffect::BehaviorRateEffect(const EffectInfo * pEffectInfo) :
	Effect(pEffectInfo)
{
	this->lpBehaviorData = 0;
	this->linitialValues = 0;
	this->lvalues = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void BehaviorRateEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	Effect::initialize(pData, pState, period, pCache);

	string name = this->pEffectInfo()->variableName();
	this->lpBehaviorData = pData->pBehaviorData(name);
	if (!this->lpBehaviorData)
	{
		throw logic_error(
			"Data for behavior variable '" + name +"' expected.");
	}
	
	this->linitialValues = this->lpBehaviorData->values(this->period());
	
	this->lvalues = pState->behaviorValues(name);
}

/**
 * Returns the number of actors of the behavior variable associated with this
 * effect.
 */
int BehaviorRateEffect::n() const
{
	return this->lpBehaviorData->n();
}


/**
 * Returns the behavior value of the given actor.
 */
int BehaviorRateEffect::value(int actor) const
{
	return this->lvalues[actor];
}


/**
 * Returns the initial behavior value of the given actor.
 */
int BehaviorRateEffect::initialValue(int actor) const
{
	return this->linitialValues[actor];
}


/**
 * Returns the overall mean of all observed values.
 */
double BehaviorRateEffect::overallCenterMean() const
{
	return this->lpBehaviorData->overallMean();
}


/**
 * Returns the behavior value of the given actor centered around the
 * overall mean of all observed values.
 */
double BehaviorRateEffect::centeredValue(int actor) const
{
	return this->lvalues[actor] - this->lpBehaviorData->overallMean();
}


/**
 * Returns if the value of the behavioral variable is missing for the given
 * actor at the specified observation.
 */
bool BehaviorRateEffect::missing(int observation, int actor) const
{
	return this->lpBehaviorData->missing(observation, actor);
}


/**
 * Returns the observed range of the respective behavior variable.
 */
double BehaviorRateEffect::range() const
{
	return this->lpBehaviorData->range();
}


/**
 * Returns the centered similarity for the given values defined as
 * 1 - |a - b| / range - similarityMean.
 */
double BehaviorRateEffect::similarity(double a, double b) const
{
	return this->lpBehaviorData->similarity(a, b);
}


/**
 * Returns the similarity mean value over all observations.
 */
double BehaviorRateEffect::similarityMean() const
{
	return this->lpBehaviorData->similarityMean();
}


/**
 * Returns the value of the variance  over all observations.
 */
double BehaviorRateEffect::variance() const
{
    return this->lpBehaviorData->variance();
}

}