/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BehaviorEffect.
 *****************************************************************************/

#include <stdexcept>

#include "ContinuousEffect.h"
#include "data/Data.h"
#include "data/ContinuousLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/ContinuousVariable.h"

namespace siena
{

/**
 * Constructor.
 */
ContinuousEffect::ContinuousEffect(const EffectInfo * pEffectInfo) :
	Effect(pEffectInfo)
{
	this->lpContinuousData = 0;
	this->lvalues = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ContinuousEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	Effect::initialize(pData, pState, period, pCache);

	string name = this->pEffectInfo()->variableName();

	this->lpContinuousData = pData->pContinuousData(name);

	if (!this->lpContinuousData)
	{
		throw logic_error(
			"Data for behavior variable '" + name +"' expected.");
	}

	this->lvalues = pState->continuousValues(name);
}


/**
 * Sets the coefficient of this effect in the Bergstrom step.
 */
void ContinuousEffect::coefficient(double value)
{
	this->lcoefficient = value;
}


/**
 * Returns the number of actors of the continuous behavior variable 
 * associated with this effect.
 */
int ContinuousEffect::n() const
{
	return this->lpContinuousData->n();
}


/**
 * Returns the continuous behavior value of the given actor.
 */
double ContinuousEffect::value(int actor) const
{
	return this->lvalues[actor];
}


/**
 * Returns the continuous behavior value of the given actor centered 
 * around the overall mean of all observed values.
 */
double ContinuousEffect::centeredValue(int actor) const
{
	return this->lvalues[actor] - this->lpContinuousData->overallMean();
}


/**
 * Returns if the value of the continuous behavior variable is missing 
 * for the given actor at the specified observation.
 */
bool ContinuousEffect::missing(int observation, int actor) const
{
	return this->lpContinuousData->missing(observation, actor);
}


/**
 * Returns the observed range of the respective continuous behavior 
 * variable.
 */
double ContinuousEffect::range() const
{
	return this->lpContinuousData->range();
}


/**
 * Returns the centered similarity for the given values defined as
 * 1 - |a - b| / range - similarityMean.
 */
double ContinuousEffect::similarity(double a, double b) const
{
	return this->lpContinuousData->similarity(a, b);
}


/**
 * Returns the similarity mean value over all observations.
 */
double ContinuousEffect::similarityMean() const
{
	return this->lpContinuousData->similarityMean();
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given values of
 * the behavior variable.
 */
double ContinuousEffect::evaluationStatistic(double * currentValues)
{
	double statistic = 0;
	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		this->preprocessEgo(i);
		if (!this->missing(this->period(), i) &&
			!this->missing(this->period() + 1, i))
		{
			statistic += this->egoStatistic(i, currentValues);
		}
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double ContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	throw runtime_error("egoStatistic not implemented for " +
		this->pEffectInfo()->effectName());
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to the given values of
 * the behavior variable.
 *//*
double BehaviorEffect::endowmentStatistic(const int * difference,
	double * currentValues)
{
	double statistic = 0;
	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		this->preprocessEgo(i);
		if (!this->missing(this->period(), i))
		{
			statistic += this->egoEndowmentStatistic(i, difference,
				currentValues);
		}
	}

	return statistic;
}*/
/**
 * Returns the statistic corresponding the given ego as part of
 * the endowment function with respect to an initial behavior
 * variable and the current state.
 */
double ContinuousEffect::egoEndowmentStatistic(int i, const int * difference,
	double *currentValues)
{
	throw runtime_error("egoEndowmentStatistic not implemented for " +
		this->pEffectInfo()->effectName());
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the creation function.
 * @param[in] difference an array of differences per each actor where the
 * current value is subtracted from the initial value. Thus positive
 * differences indicate a decrease of actors' behavior, while negative values
 * indicate an increase of actors' behavior.
 * @param[in] currentValues the current state of the behavior variable
 *//*
double BehaviorEffect::creationStatistic(int * difference,
	double *currentValues)
{
	// Here we use a trick. The creation statistics are very similar to the
	// endowmnent statistics, but instead of summing over all actors with
	// decreasing values, we must now sum over all actors with increasing
	// values. So we just reverse the differences and call the endowment
	// statistic.

	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		difference[i]  = -difference[i];
	}

	double statistic = this->endowmentStatistic(difference, currentValues);

	for (int i = 0; i < n; i++)
	{
		difference[i]  = -difference[i];
	}

	return -statistic;
}
*/

/**
 * Does the necessary preprocessing work for calculating the probabilities
 * for a specific ego. This method must be invoked before
 * calling BehaviorEffect::calculateChangeContribution(...).
 */
void ContinuousEffect::preprocessEgo(int ego)
{
	this->lego = ego;
}
}
