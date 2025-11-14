/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterCcEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageAlterCcEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageAlterCcEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AverageAlterCcEffect::AverageAlterCcEffect(
	const EffectInfo * pEffectInfo, bool divide) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
}

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] reciprocal indicates if only reciprocal ties have to be
 * considered
 * @param simulatedState If `true` the value() function uses the simulated
 *        state, if any or the value at the end of the period.
 */
 AverageAlterCcEffect::AverageAlterCcEffect(
	const EffectInfo * pEffectInfo, bool divide, const bool simulatedState) :
		NetworkDependentBehaviorEffect(pEffectInfo, simulatedState)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AverageAlterCcEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();
	double currentMean = 0; // contemporaneous mean, inefficiently calculated each time here
	
	contribution = totalAlterValue(actor); // this is sum of grand-mean centered alter values
	                                       // and zero if there are not alters;
										   // see NetworkDependentBehaviorEffect.cpp
	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = v_i * avg(v_j) over all neighbors j of i.
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// This is d * avg(v_j), the average being taken over all neighbors
		// of i. This is what is calculated below.
		// if (not divide), instead of avg the total is used.
		
		// calculate contemporaneous mean to re-center
		// the grand-mean centered value
		for (int i = 0; i < this->n(); i++)
		{
			currentMean += this->value(i);
		}
		currentMean /= this->n();

		// re-center totalAlterValue
		contribution += pNetwork->outDegree(actor) * (overallCenterMean() - currentMean);
		// weight by step
		contribution *= difference;

		if (this->ldivide)
		{
			contribution /= pNetwork->outDegree(actor);
		}
	}

	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AverageAlterCcEffect::egoStatistic(int i, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;
	double currentDelta = 0; // correction for contemporaneous mean, inefficiently calculated each time here

	for (IncidentTieIterator iter = pNetwork->outTies(i);
		 iter.valid();
		 iter.next())
	{
		statistic += currentValues[iter.actor()];
		neighborCount++;
	}

	if (neighborCount > 0)
	{
		
		// calculate difference between contemporaneous mean and grand mean
		for (int j = 0; j < this->n(); j++)
		{
			currentDelta += currentValues[j]; // note these are centered at grand mean
		}
		currentDelta /= this->n();

		// re-center current statistic value to current mean
		statistic -= neighborCount * currentDelta;
		// multiply with re-centered ego value
		statistic *= (currentValues[i] - currentDelta);
		if (this->ldivide)
		{
			statistic /= neighborCount;
		}
	}

	return statistic;
}

}
