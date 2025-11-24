/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreeWeightedAverageGroupEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IndegreeWeightedAverageGroupEffect class.
 *****************************************************************************/

//#include <R_ext/Print.h>

/**
 * Not all includes might be necessary.
*/

#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "IndegreeWeightedAverageGroupEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
IndegreeWeightedAverageGroupEffect::IndegreeWeightedAverageGroupEffect(const EffectInfo * pEffectInfo, bool divide) :
	NetworkDependentBehaviorEffect(pEffectInfo)
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
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void IndegreeWeightedAverageGroupEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkDependentBehaviorEffect::initialize(pData, pState, period, pCache);
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double IndegreeWeightedAverageGroupEffect::calculateChangeContribution(int actor,
	int difference)
{
	double statistic = 0;
    int weightedN = 0;
/* 	const Network * pNetwork = this->pNetwork();*/
	for (int i = 0; i < this->n(); i++)
	{
		statistic += this->centeredValue(i) * this->pNetwork()->inDegree(i);
		if (this->ldivide)
		{
	        weightedN += this->pNetwork()->inDegree(i);
		}
	}
	statistic += this->centeredValue(actor) * this->pNetwork()->inDegree(actor) + difference;
	if (this->ldivide)
	{
	    weightedN += this->pNetwork()->inDegree(actor);
	    statistic /= weightedN;
	}
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
double IndegreeWeightedAverageGroupEffect::egoStatistic(int ego, double * currentValues)
{
	double thesum = 0;
    int weightsum = 0;
 	for (int i = 0; i < this->n(); i++)
	{
		thesum += currentValues[i] * this->pNetwork()->inDegree(i);
		if (this->ldivide)
		{
        weightsum += this->pNetwork()->inDegree(i);
		}
	}
	if (this->ldivide)
	{
	    thesum /= weightsum;    
	}
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
double IndegreeWeightedAverageGroupEffect::egoEndowmentStatistic(int ego,
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
