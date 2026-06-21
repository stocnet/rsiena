/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DegreeWeightedAverageGroupEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DegreeWeightedAverageGroupEffect class.
 *****************************************************************************/

//#include <R_ext/Print.h>

/**
 * Not all includes might be necessary.
*/

#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "DegreeWeightedAverageGroupEffect.h"
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
DegreeWeightedAverageGroupEffect::DegreeWeightedAverageGroupEffect(const EffectInfo * pEffectInfo, 
		bool divide, bool outdegree, bool nc, bool ego) :
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
	this->loutdegree = outdegree;
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
void DegreeWeightedAverageGroupEffect::initialize(const Data * pData,
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
double DegreeWeightedAverageGroupEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
    int weightedN = 0;
	double multiplier = 1.0;

/* 	const Network * pNetwork = this->pNetwork();*/
	for (int i = 0; i < this->n(); i++)
	{
		if (i == actor && !this->lego)
		{
			continue; // Pure normative context: skip ego!
		}
		int degree = 0;
		degree = this->loutdegree ? this->pNetwork()->outDegree(i) :
										   this->pNetwork()->inDegree(i);
		contribution += resolvedValue(i, this->lnc) * degree;
		if (this->ldivide || !this->lcenterMean)
		{
	        weightedN += degree;
		}
	}
	int actorDegree = 0;
	actorDegree = this->loutdegree ? this->pNetwork()->outDegree(actor) :
										  this->pNetwork()->inDegree(actor);
	if (this->lego)
	{
		contribution += (resolvedValue(actor, this->lnc) + difference) * actorDegree;
	}
	if (this->ldivide)
	{
	    contribution /= weightedN;
	}
	else if (!this->lcenterMean)
	{
		multiplier = weightedN;
	}
	// recentering is not meaningful for the uncentered branch
	if (!this->lcenterMean && !this->lnc)
	{
		contribution += multiplier * (this->overallCenterMean() - this->lcenteringValue);
	}
	contribution *= difference;
	return contribution; 
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double DegreeWeightedAverageGroupEffect::egoStatistic(int ego, 
	double * currentValues)
{
	double statistic = 0;
    int weightsum = 0;
	double multiplier = 1.0;
 	for (int i = 0; i < this->n(); i++)
	{
		if (!this->lego && i == ego)
		{
			continue;
		}
		double value = this->lnc ? currentValues[i] + this->overallCenterMean() : 
								   currentValues[i];
		double degree = this->loutdegree ? this->pNetwork()->outDegree(i) : 
										   this->pNetwork()->inDegree(i);
		statistic += value * degree;
		if (this->ldivide || !this->lcenterMean)
		{
			weightsum += degree;
		}
	}
	if (this->ldivide)
	{
	    statistic /= weightsum;    
	}
	else if (!this->lcenterMean)
	{
		multiplier = weightsum;
	}
	if (!this->lnc && !this->lcenterMean) // recentering makes no cense if non-centered values are used	
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
double DegreeWeightedAverageGroupEffect::egoEndowmentStatistic(int ego,
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
