/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAverageEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AltersCovariateAverageEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "AltersCovariateAverageEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
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
AltersCovariateAverageEffect::AltersCovariateAverageEffect(
	const EffectInfo * pEffectInfo, bool divide) :
	CovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
}



/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AltersCovariateAverageEffect::calculateChangeContribution(int actor,
	int difference)
{
	double statistic = 0;
	if (this->ldivide)
	{
		statistic = difference * this->averageAlterValue(actor);
	}
	else
	{
		statistic = difference * this->totalAlterValue(actor);
	}
	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AltersCovariateAverageEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;

	if (!this->missingDummy(ego))
	{
		if (this->ldivide)
		{
			statistic = currentValues[ego] * this->averageAlterValue(ego);
		}
		else
		{
			statistic = currentValues[ego] * this->totalAlterValue(ego);
		}
	}
	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AltersCovariateAverageEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if (difference[ego] > 0 && !this->missingDummy(ego))
	{
		if (this->ldivide)
		{
			statistic -= difference[ego] * this->averageAlterValue(ego);
		}
		else
		{
			statistic -= difference[ego] * this->totalAlterValue(ego);
		}
	}
	return statistic;
}


}
