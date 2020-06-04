/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InterceptEffect.cpp
 *
 * Description: This file contains the implementation of the
 * InterceptEffect class.
 *****************************************************************************/

#include "InterceptEffect.h"
#include "model/variables/ContinuousVariable.h"
#include "data/ContinuousLongitudinalData.h"

namespace siena
{

/**
 * Constructor.
 */
InterceptEffect::InterceptEffect(const EffectInfo * pEffectInfo) :
	ContinuousEffect(pEffectInfo)
{
}


/**
 * Returns how much this effect contributes to the change in the
 * continuous behavior.
 */
double InterceptEffect::calculateChangeContribution(int actor)
{
	return 1;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double InterceptEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego];
}


/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
//double InterceptEffect::egoEndowmentStatistic(int ego,
//	const int * difference,
//	double * currentValues)
//{
//	double statistic = 0;
//
//	if (difference[ego] > 0)
//	{
//		//	statistic = currentValues[ego] ;
//		statistic -= difference[ego];
//	}
//	return statistic;
//}

}
