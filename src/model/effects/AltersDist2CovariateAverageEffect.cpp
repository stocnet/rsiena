/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersDist2CovariateAverageEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AltersDist2CovariateAverageEffect class.
 * It is like AverageAlterDist2Effect,
 * but now an extension of CovariateAndNetworkBehaviorEffect.
 *****************************************************************************/

#include "AltersDist2CovariateAverageEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
AltersDist2CovariateAverageEffect::AltersDist2CovariateAverageEffect(
	const EffectInfo * pEffectInfo, bool divide1, bool divide2) :
	CovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->ldivide1 = divide1;
	// Indicates whether there will be division by the outdegree of ego
	this->ldivide2 = divide2;
	// Indicates whether there will be division by the outdegree of alter
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AltersDist2CovariateAverageEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = v_i * avg(v_j) over all neighbors j of i,
		// where v_j is the average covariate of j's neighbors,
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
			double alterValue = this->totalAlterValue(iter.actor());
			int tieValue =  this->pNetwork()->tieValue(iter.actor(), actor);
			if (tieValue == 1)
			{
				alterValue -= this->centeredValue(actor);
			}
			if (this->ldivide2)
			{
				if ((pNetwork->outDegree(iter.actor()) - tieValue)> 0)
				{
					alterValue /= (pNetwork->outDegree(iter.actor()) - tieValue);
				}
				else
				{
					alterValue = this->covariateMean();
				}
			}
			sumAlterValue += alterValue;
		}
		contribution = difference * sumAlterValue;
		if (this->ldivide1)
		{
			contribution /= pNetwork->outDegree(actor);
		}
	}
	else
	{
		if (this->ldivide1)
		{
			contribution = this->covariateMean();
		}
	}

	return contribution;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AltersDist2CovariateAverageEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();
		double alterXValue = 0;
		int tieToEgo = 0;
		for (IncidentTieIterator iteri = pNetwork->outTies(j);
			iteri.valid();
			iteri.next())
		{
			if (ego != iteri.actor())
			{
				alterXValue += this->covariateValue(iteri.actor());
			}
			else
			{
				tieToEgo = 1;
			}
		}
// tieToEgo =  this->pNetwork()->tieValue(iter.actor(), ego);
		if (this->ldivide2)
		{
			if (pNetwork->outDegree(j) > tieToEgo)
			{
				alterXValue /= (pNetwork->outDegree(j) - tieToEgo);
			}
			else
			{
				alterXValue = this->covariateMean();
			}
		}
		statistic += alterXValue;
		neighborCount++;
	}

	if (neighborCount > 0)
	{
		statistic *= currentValues[ego];
		if (this->ldivide1)
		{
			statistic /= neighborCount;
		}
	}
	else
	{
		if (this->ldivide1)
		{
			statistic = this->covariateMean();
		}
	}

	return statistic;
}

}
