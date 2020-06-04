/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimilarityWEffect.cpp
 *
 * Description: This file contains the implementation of the
 * SimilarityWEffect class.
 *****************************************************************************/
#include <cstdlib>
#include "SimilarityWEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] average indicates if one of the average effects is required
 * @param[in] alterPopularity indicates if the similarity scores have to
 * be multiplied by the in-degrees of alters
 */
SimilarityWEffect::SimilarityWEffect(
		const EffectInfo * pEffectInfo,
		bool average,
		bool alterPopularity,
		bool egoPopularity) :
	DyadicCovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->laverage = average;
	this->lalterPopularity = alterPopularity;
	this->legoPopularity = egoPopularity;
	this->lpar2 = (pEffectInfo->internalEffectParameter() >= 2);
	// specifies type of denominator
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double SimilarityWEffect::calculateChangeContribution(int actor,
		int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the average similarity W effect:
		// s_i(x) = avg(w(i,j)*(sim(v_i, v_j) - centeringConstant)) over all neighbors
		// j of i.
		// sim(v_i, v_j) = 1.0 - |v_i - v_j| / observedRange
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// To this end, we can disregard the centering constant and
		// compute the average change in similarity, namely,
		// avg(w(i,j)*(sim(v_i + d, v_j) - sim(v_i, v_j))) =
		// avg(w(i,j)*(1 - |v_i+d-v_j|/range - 1 + |v_i-v_j|/range)) =
		// avg(w(i,j)*(|v_i-v_j| - |v_i+d-v_j|)) / range,
		// the average being taken over all neighbors of i.
		// The reasoning for avg. similarity x popularity alter W effect is
		// similar.
		// This is what is calculated below.

		int oldValue = this->value(actor);
		int newValue = oldValue + difference;
		int totalCount = 0;
		double totalChange = 0;
		double totalWeightValue = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();
			int alterValue = this->value(j);
			double dycova = this->dycoValue(actor, j);
			double change = dycova *
				(std::abs(oldValue - alterValue) - std::abs(newValue - alterValue));
			if (dycova != 0.0) totalCount++;
			if (this->laverage)
			{
				if (lpar2)
				{
					totalWeightValue += dycova;
				}
				else
				{
					totalWeightValue += 1;
				}
			}

			if (this->lalterPopularity)
			{
				change *= pNetwork->inDegree(j);
			}

			totalChange += change;
		}

		contribution = totalChange / this->range();

		if ((this->laverage) && (totalCount > 0))
		{
			contribution /= totalWeightValue;
		}

		if (this->legoPopularity)
		{
			contribution *= pNetwork->inDegree(actor);
		}
	}
	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double SimilarityWEffect::egoStatistic(int ego,
	double * currentValues)
{
	const Network * pNetwork = this->pNetwork();

	double statistic = 0;
	double totalWeightValue = 0;
	int totalCount = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j))
		{
			double dycova = this->dycoValue(ego, j);
			double tieStatistic = dycova *
				this->similarity(currentValues[ego], currentValues[j]);
			if (this->lalterPopularity)
			{
				tieStatistic *= pNetwork->inDegree(j);
			}
			statistic += tieStatistic;
			if (dycova != 0.0) totalCount++;
			if (lpar2)
			{
				totalWeightValue += dycova;
			}
			else
			{
				totalWeightValue += 1;
			}
		}
	}
		if ((this->laverage) && (totalCount > 0))
	{
		statistic /= totalWeightValue;
	}

	if (this->legoPopularity)
	{
		statistic *= pNetwork->inDegree(ego);
	}

	return statistic;
}

}
