/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterInDist2Effect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageAlterInDist2Effect class.
 *****************************************************************************/

#include <cmath>
#include "AverageAlterInDist2Effect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "NetworkDependentBehaviorEffect.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AverageAlterInDist2Effect::AverageAlterInDist2Effect(
	const EffectInfo * pEffectInfo, bool divide1, bool divide2) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldivide1 = divide1;
	// Indicates whether there will be division by the outdegree of ego
	this->ldivide2 = divide2;
	// Indicates whether there will be division by the indegree of alter
}

/**
 * Deallocates this effect object;
 */
AverageAlterInDist2Effect::~AverageAlterInDist2Effect()
{
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 * It is assumed that preprocessEgo(ego) has been called before.

 */
double AverageAlterInDist2Effect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = v_i * avg(v_j) over all neighbors j of i,
		// where v_j is the average behavior of j's in-neighbors,
		// excluding i.
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// This is d * avg(v_j) and is calculated below.
		// if (not divide1) or (not divide2),
		// instead of "avg" or "average" the total is used.

		double sumAlterValue = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			double alterValue = this->totalInAlterValue(iter.actor());
			int tieValue =  this->pNetwork()->tieValue(actor, iter.actor());
			if (tieValue == 1)
			{
				alterValue -= this->centeredValue(actor);
			}
			if (((pNetwork->inDegree(iter.actor()) - tieValue)> 0) & (this->ldivide2))
			{
				alterValue /= (pNetwork->inDegree(iter.actor()) - tieValue);
			}
			sumAlterValue += alterValue;
		}
		contribution = difference * sumAlterValue;
		if (this->ldivide1)
		{
			contribution /= pNetwork->outDegree(actor);
		}
	}
	return contribution;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * currentValues given for the behavior variable.
 */
double AverageAlterInDist2Effect::egoStatistic(int i, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(i);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();
		double alterValue = 0;
		for (IncidentTieIterator iteri = pNetwork->inTies(j);
			iteri.valid();
			iteri.next())
		{
			if (i != iteri.actor())
			{
				alterValue += currentValues[iteri.actor()];
			}
		}
// tieFromi =  this->pNetwork()->tieValue(i, iter.actor());
		if ((pNetwork->inDegree(j) > 1) & (this->ldivide2))
		{
			alterValue /= (pNetwork->inDegree(j) - 1);
				// there always is a tie i -> iteri.actor()
		}
		statistic += alterValue;
		neighborCount++;
	}

	if (neighborCount > 0)
	{
		statistic *= currentValues[i];
		if (this->ldivide1)
		{
			statistic /= neighborCount;
		}
	}

	return statistic;
}

}
