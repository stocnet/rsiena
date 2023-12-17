/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegWeightAverageEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IndegWeightAverageEffect class.
 *****************************************************************************/

//#include <R_ext/Print.h>

#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "IndegWeightAverageEffect.h"
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
IndegWeightAverageEffect::IndegWeightAverageEffect(const EffectInfo * pEffectInfo) :
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
}

/**
 * Might not be correct to do
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void IndegWeightAverageEffect::initialize(const Data * pData,
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
double IndegWeightAverageEffect::calculateChangeContribution(int actor,
	int difference)
{
	double statistic = 0;
    int weightedN = 0;
/* 	const Network * pNetwork = this->pNetwork();*/
	for (int i = 0; i < this->n(); i++)
	{
		statistic += this->centeredValue(i) * this->pNetwork()->inDegree(i);
        weightedN += this->pNetwork()->inDegree(i);
	}
	statistic += this->centeredValue(actor) * this->pNetwork()->inDegree(actor) + difference;
    weightedN += this->pNetwork()->inDegree(actor);
    statistic /= weightedN;
/**
 * the centering might not work as intended anymore; use only for p = 0 for now
*/
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
double IndegWeightAverageEffect::egoStatistic(int ego, double * currentValues)
{
	double thesum = 0;
    int weightsum = 0;
 	for (int i = 0; i < this->n(); i++)
	{
		thesum += currentValues[i] * this->pNetwork()->inDegree(i);
        weightsum += this->pNetwork()->inDegree(i);
	}
	/* thesum /= this->n(); */
    thesum /= weightsum;    
   /*
    * the centering might not work as intended anymore; use only for p = 0 for now
    */
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
double IndegWeightAverageEffect::egoEndowmentStatistic(int ego,
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
