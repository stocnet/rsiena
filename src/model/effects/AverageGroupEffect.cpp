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
AverageGroupEffect::AverageGroupEffect(const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo)
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
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AverageGroupEffect::calculateChangeContribution(int actor,
	int difference)
{
	double statistic = 0;

	for (int i = 0; i < this->n(); i++)
	{
		statistic += this->centeredValue(i);
	}
	statistic += this->centeredValue(actor) + difference;
	statistic /= this->n();
//Rprintf("calculateSumPlus %f ", thesum);
//Rprintf(" and %f \n", thesum);
	if (!this->lcenterMean)
	{
		statistic += (this->overallCenterMean() - this->lcenteringValue);
	}
	return difference * statistic;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AverageGroupEffect::egoStatistic(int ego, double * currentValues)
{
	double thesum = 0;
	for (int i = 0; i < this->n(); i++)
	{
		thesum += currentValues[i];
	}
	thesum /= this->n();
	if (!this->lcenterMean)
	{
		thesum += (this->overallCenterMean() - this->lcenteringValue);
	}
	return thesum * currentValues[ego];
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
