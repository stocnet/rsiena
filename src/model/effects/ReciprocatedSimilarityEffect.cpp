/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocatedSimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ReciprocatedSimilarityEffect class.
 *****************************************************************************/

#include <cstdlib>
#include <stdexcept>
#include <string>
#include "ReciprocatedSimilarityEffect.h"
#include "network/OneModeNetwork.h"
#include "network/CommonNeighborIterator.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] average indicates if one of the average effects is required
 * @param[in] alterPopularity indicates if the similarity scores have to
 * be multiplied by the in-degrees of alters
 */
ReciprocatedSimilarityEffect::ReciprocatedSimilarityEffect(
	const EffectInfo * pEffectInfo,
	bool average,
	bool alterPopularity) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->laverage = average;
	this->lalterPopularity = alterPopularity;
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double ReciprocatedSimilarityEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(string("One-mode network expected in ") +
			"ReciprocatedSimilarityEffect");
	}

	if (pNetwork->reciprocalDegree(actor) > 0)
	{
		// The formula for the average similarity x reciprocity effect:
		// s_i(x) = avg(sim(v_i, v_j) - centeringConstant) over all
		// mutual neighbors j of i.
		// sim(v_i, v_j) = 1.0 - |v_i - v_j| / observedRange
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// To this end, we can disregard the centering constant and
		// compute the average change in similarity, namely,
		// avg(sim(v_i + d, v_j) - sim(v_i, v_j)) =
		// avg(1 - |v_i+d-v_j|/range - 1 + |v_i-v_j|/range) =
		// avg(|v_i-v_j| - |v_i+d-v_j|) / range,
		// the average being taken over all mutual neighbors of i.
		// The reasoning for other effects is similar.
		// This is what is calculated below.

		int oldValue = this->value(actor);
		int newValue = oldValue + difference;
		int totalChange = 0;

		for (CommonNeighborIterator iter = pNetwork->reciprocatedTies(actor);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();
			int alterValue = this->value(j);
			int change =
				abs(oldValue - alterValue) - abs(newValue - alterValue);

			if (this->lalterPopularity)
			{
				change *= pNetwork->inDegree(j);
			}

			totalChange += change;
		}

		contribution = ((double) totalChange) / this->range();

		if (this->laverage)
		{
			contribution /= pNetwork->reciprocalDegree(actor);
		}
	}

	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double ReciprocatedSimilarityEffect::egoStatistic(int ego,
	double * currentValues)
{
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(string("One-mode network expected in ") +
			"ReciprocatedSimilarityEffect");
	}

	double statistic = 0;
	int neighborCount = 0;

	for (CommonNeighborIterator iter = pNetwork->reciprocatedTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j))
		{
			double tieStatistic =
				this->similarity(currentValues[ego], currentValues[j]);

			if (this->lalterPopularity)
			{
				tieStatistic *= pNetwork->inDegree(j);
			}

			statistic += tieStatistic;
			neighborCount++;
		}
	}

	if (this->laverage && neighborCount > 0)
	{
		statistic /= neighborCount;
	}

	return statistic;
}

}
