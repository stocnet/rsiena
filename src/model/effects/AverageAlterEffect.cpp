/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageAlterEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageAlterEffect.h"
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
AverageAlterEffect::AverageAlterEffect(
	const EffectInfo * pEffectInfo, bool divide, bool alterPopularity) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
	this->lalterPopularity = alterPopularity;
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AverageAlterEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	contribution = totalAlterValue(actor);
	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = v_i * avg(v_j) over all neighbors j of i.
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// This is d * avg(v_j), the average being taken over all neighbors
		// of i. This is what is calculated below.
		// if (not divide), instead of avg the total is used.
		if (this->lalterPopularity)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(actor);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				contribution += (this->centeredValue(j)) * (pNetwork->inDegree(j));
			}
		}
		else
		{
		contribution = difference * totalAlterValue(actor);
		}
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
double AverageAlterEffect::egoStatistic(int i, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(i);
		 iter.valid();
		 iter.next())
	{
		if (this->lalterPopularity)
		{
			int j = iter.actor();
			statistic += (currentValues[j]) * (pNetwork->inDegree(j));
		}
		else
		{
			statistic += currentValues[iter.actor()];
		}
		neighborCount++;
	}

	if (neighborCount > 0)
	{
		statistic *= currentValues[i];
		if (this->ldivide)
		{
			statistic /= neighborCount;
		}
	}

	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AverageAlterEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	const Network * pNetwork = this->pNetwork();

	if (difference[ego] > 0)
	{
		if (pNetwork->outDegree(ego) > 0)
		{
			double thisStatistic = 0;
			double previousStatistic = 0;

			for (IncidentTieIterator iter = pNetwork->outTies(ego);
				 iter.valid();
				 iter.next())
			{
				double alterValue = currentValues[iter.actor()];
				double alterPreviousValue = currentValues[iter.actor()]
													+ difference[iter.actor()]; // this addition was absent until 1.1-292
				thisStatistic += alterValue;
				previousStatistic += alterPreviousValue;
//				thisStatistic += iter.value() * alterValue; these lines also until 1.1-291; but iter.value is always 1?
//				previousStatistic += iter.value() * alterPreviousValue;
			}

			thisStatistic *= currentValues[ego];
			previousStatistic *= (currentValues[ego] + difference[ego]);
			statistic = thisStatistic - previousStatistic;
			if (this->ldivide)
			{
				statistic /= pNetwork->outDegree(ego);
			}
		}
	}

	return statistic;
}

}
