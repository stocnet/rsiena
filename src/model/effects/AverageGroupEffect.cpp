/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageGroupEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageGroupEffect class.
 *****************************************************************************/

//#include <R_ext/Print.h>

#include "AverageGroupEffect.h"
#include "model/variables/BehaviorVariable.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/EffectInfo.h"

namespace siena
{

/**
 * Constructor.
 */
AverageGroupEffect::AverageGroupEffect(const EffectInfo * pEffectInfo, bool divide,
	bool nc, bool ego) : BehaviorEffect(pEffectInfo)
{
	this->lcenterMean = (pEffectInfo->internalEffectParameter() <= 0.5);
	if (!this->lcenterMean)
	{
		this->lcenteringValue = pEffectInfo->internalEffectParameter();
	}
	else
	{
		this->lcenteringValue = 0.0;
	}
	this->ldivide = divide;
	// Indicates whether there will be division by the number of actors.
	this->lnc = nc;
	this->lego = ego;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void AverageGroupEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	BehaviorEffect::initialize(pData, pState, period, pCache);
	/// calculate groupsize that is constant over egos
	this->lgroupSize = !this->lego ? this->n() - 1 : this->n();
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AverageGroupEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	double multiplier = 1.0;

	for (int i = 0; i < this->n(); i++)
	{
		if (!this->lego && i == actor)
		{
			continue; // Pure normative context: skip ego!
		}
		contribution += this->resolvedValue(i, this->lnc);
	}
	if (this->lego)
	{
		contribution += this->resolvedValue(actor, this->lnc) + difference;
	}
	if (this->ldivide) 
	{
		contribution /= lgroupSize;
	}
	else if (!this->lcenterMean)
	{
		multiplier = lgroupSize;
	}
//Rprintf("calculateSumPlus %f ", thesum);
//Rprintf(" and %f \n", thesum);
	// recentering only makes sense for the centered branch
	if (!this->lcenterMean && !this->lnc)
	{
		// multiply by n for the sum to scale correctly for total change
		contribution += multiplier * (this->overallCenterMean() - this->lcenteringValue);
	}
	contribution *= difference;
	return contribution;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AverageGroupEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	double multiplier = 1.0;
    for (int i = 0; i < this->n(); i++)
    {
        if (!this->lego && i == ego)
		{
			continue;
		}
        statistic += this->lnc ? currentValues[i] + this->overallCenterMean() : 
							     currentValues[i];
    }
	if (this->ldivide)
	{
		statistic /= this->lgroupSize;
	}
	else if (!this->lcenterMean)
	{
		multiplier = this->lgroupSize;
	}
	// recentering only makes sense for the centered branch
	if (!this->lcenterMean && !this->lnc)
	{
		statistic += multiplier * (this->overallCenterMean() - this->lcenteringValue);
	}
	statistic *= this->lnc ? currentValues[ego] + this->overallCenterMean()  : 
							 currentValues[ego];
	return statistic;
}


/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AverageGroupEffect::egoEndowmentStatistic(int ego,
 	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if (difference[ego] > 0)
	{
		statistic = difference[ego] * this->egoStatistic(ego, currentValues);
	}

	return statistic;
}

}
