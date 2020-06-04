/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorInteractionEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BehaviorInteractionEffect.
 *****************************************************************************/

#include "BehaviorInteractionEffect.h"

namespace siena
{

/**
 * Constructs a new interaction effect between the given effects.
 * The parameter pEffect3 should be 0 for two-way interactions.
 * This effect takes the ownership of the given effects, which mean
 * that the given effects will be destroyed as soon as this
 * effect is destroyed.
 */
BehaviorInteractionEffect::BehaviorInteractionEffect(
	const EffectInfo * pEffectInfo,
	BehaviorEffect * pEffect1,
	BehaviorEffect * pEffect2,
	BehaviorEffect * pEffect3) : BehaviorEffect(pEffectInfo)
{
	this->lpEffect1 = pEffect1;
	this->lpEffect2 = pEffect2;
	this->lpEffect3 = pEffect3;
}


/**
 * Deallocates this effects.
 */
BehaviorInteractionEffect::~BehaviorInteractionEffect()
{
	delete this->lpEffect1;
	delete this->lpEffect2;
	delete this->lpEffect3;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void BehaviorInteractionEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	BehaviorEffect::initialize(pData, pState, period, pCache);
	this->lpEffect1->initialize(pData, pState, period, pCache);
	this->lpEffect2->initialize(pData, pState, period, pCache);

	if (this->lpEffect3)
	{
		this->lpEffect3->initialize(pData, pState, period, pCache);
	}
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling BehaviorEffect::calculateChangeContribution(...).
 */
void BehaviorInteractionEffect::preprocessEgo(int ego)
{
	BehaviorEffect::preprocessEgo(ego);

	this->lpEffect1->preprocessEgo(ego);
	this->lpEffect2->preprocessEgo(ego);

	if (this->lpEffect3)
	{
		this->lpEffect3->preprocessEgo(ego);
	}
}


/**
 * Calculates the change in the statistic corresponding to this effecf if
 * the given actor were to change his behavior by the given amount.
 */
double BehaviorInteractionEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution =
		this->lpEffect1->calculateChangeContribution(actor, difference) *
		this->lpEffect2->calculateChangeContribution(actor, difference);

	contribution /= difference;

	if (this->lpEffect3)
	{
		contribution *= this->lpEffect3->calculateChangeContribution(actor,
			difference) / difference;
	}

	return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double BehaviorInteractionEffect::egoStatistic(int i,
	double * currentValues)
{
	double statistic;

	if (currentValues[i] != 0)
	{
		statistic = this->lpEffect1->egoStatistic(i, currentValues) *
			this->lpEffect2->egoStatistic(i, currentValues);

		statistic /= currentValues[i];

		if (this->lpEffect3)
		{
			statistic *= this->lpEffect3->egoStatistic(i, currentValues)
				/ currentValues[i];
		}
	}
	else
	{
		statistic = 0;
	}

	return statistic;
}
/**
 * Calculates the statistic corresponding to the given ego with respect to the
 * initial values of the behavior variable and the current values.
 */
double BehaviorInteractionEffect::egoEndowmentStatistic(int i,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if (difference[i] > 0)
	{
		statistic = this->lpEffect1->egoEndowmentStatistic(i,
			difference, currentValues) *
			this->lpEffect2->egoEndowmentStatistic(i, difference,
				currentValues);

		//	statistic /= currentValues[i];
		statistic /=  - difference[i];

		if (this->lpEffect3)
		{
			// statistic *= this->lpEffect3->egoEndowmentStatistic(i, difference,
			// 	currentValues) / currentValues[i];
			statistic *= this->lpEffect3->egoEndowmentStatistic(i, difference,
				currentValues) / - difference[i];
		}
	}
	else
	{
		statistic = 0;
	}

	return statistic;
}

}
