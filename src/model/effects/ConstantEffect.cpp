/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ConstantEffect class.
 *****************************************************************************/

#include <cmath>
#include "ConstantEffect.h"
#include "model/variables/BehaviorVariable.h"
#include "data/BehaviorLongitudinalData.h"

namespace siena
{

/**
 * Constructor.
 */
ConstantEffect::ConstantEffect(const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo)
{
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ConstantEffect::initialize(const Data * pData,
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
double ConstantEffect::calculateChangeContribution(int actor,
	int difference)
{	
	int statistic = -1;
	if (difference == 0)
	{
		statistic = 0;
	}
	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double ConstantEffect::egoStatistic(int ego, double * currentValues)
{
	return(- fabs(double (currentValues[ego] + 
					this->overallCenterMean() - this->initialValue(ego))));
}

}
