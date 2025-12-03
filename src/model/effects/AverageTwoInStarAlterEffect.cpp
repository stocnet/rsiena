/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageTwoInStarAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageTwoInStarAlterEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageTwoInStarAlterEffect.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
// #include "network/CommonNeighborIterator.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "NetworkDependentBehaviorEffect.h"
// #include "NetworkEffect.h"
#include "model/tables/ConfigurationTable.h"
#include "model/EffectInfo.h"

// #include "model/tables/NetworkCache.h"
// #include "model/tables/EgocentricConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AverageTwoInStarAlterEffect::AverageTwoInStarAlterEffect(
	const EffectInfo * pEffectInfo, bool divide1, bool divide2) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldivide1 = divide1;
	// Indicates whether there will be division by the outdegree of ego
	this->ldivide2 = divide2;
	// Indicates whether there will be division by the indegree of alter
	// this->lcentered = centered;
}

/**
 * Deallocates this effect object;
 */
AverageTwoInStarAlterEffect::~AverageTwoInStarAlterEffect()
{
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 * It is assumed that preprocessEgo(ego) has been called before.

 */
double AverageTwoInStarAlterEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();
	
	if (pNetwork->outDegree(actor) > 0) 
	{
		// The formula for the effect:
		// s_i(x) = v_i * sum(v_h) over all two-in-star alters h of i,
		// weighted by the number of common two-in-stars divided
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).

		double sumAlterValue = 0;
		// double denom = 0;
		for (int h = 0; h < this->n(); h++) //inefficient?
		{
			if (h != actor)
			{
				// int instars = pNetwork->inTwoStarCount(actor, j);
				int instars = this->pInStarTable()->get(h);
				double alterValue = this->value(h) * instars;
				int tieValue =  this->pNetwork()->tieValue(h, actor);
				if (((pNetwork->outDegree(h) - tieValue)> 0) && (this->ldivide2))
				{
					alterValue /= (pNetwork->outDegree(h) - tieValue);
				}
				sumAlterValue += alterValue;
			}
		}
		contribution = difference * sumAlterValue;
		// if (denom != 0)
		// {
		// 	contribution /= denom; //what happens if denom == 0 but sumaAlterValue != 0 ?
		// }
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
double AverageTwoInStarAlterEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	double alterValue = 0;

	// double denom = 0;
	for (int h = 0; h < this->n(); h++)
	{
		if (h != ego)
		{
			int instarCount = pNetwork->inTwoStarCount(ego, h);
			double alterValue = (currentValues[h] + this->overallCenterMean()) * instarCount;
			if ((pNetwork->outDegree(h) > 1) && (this->ldivide2))
			{
				alterValue /= (pNetwork->outDegree(h) - 1);
					// there always is a tie i -> iteri.actor()
			}

			// Configuration Table can not be used because preprocess ego has not been called before?
			// if (instarCount > 0)
			// {
			// 	denom += currentValues[j] + this->overallCenterMean();
			// }
			statistic += alterValue;
		}
		// tieFromi =  this->pNetwork()->tieValue(i, iter.actor());

	}
	statistic *= (currentValues[ego] + this->overallCenterMean());
	// if (denom != 0)
	// {
	// 	statistic /= denom; //what happens if denom == 0 but statistic != 0 ?
	// }
	if ((pNetwork->outDegree(ego) > 0) && (this->ldivide1))
	{
		statistic /= pNetwork->outDegree(ego);
	}
	
	return statistic;

}

}
