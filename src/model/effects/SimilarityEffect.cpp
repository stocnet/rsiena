/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * SimilarityEffect class.
 *****************************************************************************/
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "SimilarityEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] average indicates if one of the average effects is required
 * @param[in] alterPopularity indicates if the similarity scores have to
 * be multiplied by the in-degrees of alters
 */
SimilarityEffect::SimilarityEffect(
	const EffectInfo * pEffectInfo,
	bool average,
	bool alterPopularity,
	bool egoPopularity) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->laverage = average;
	this->lalterPopularity = alterPopularity;
	this->legoPopularity = egoPopularity;
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double SimilarityEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the average similarity effect:
		// s_i(x) = avg(sim(v_i, v_j) - centeringConstant) over all neighbors
		// j of i.
		// sim(v_i, v_j) = 1.0 - |v_i - v_j| / observedRange
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// To this end, we can disregard the centering constant and
		// compute the average change in similarity, namely,
		// avg(sim(v_i + d, v_j) - sim(v_i, v_j)) =
		// avg(1 - |v_i+d-v_j|/range - 1 + |v_i-v_j|/range) =
		// avg(|v_i-v_j| - |v_i+d-v_j|) / range,
		// the average being taken over all neighbors of i.
		// The reasoning for avg. similarity x popularity alter effect is
		// similar.
		// This is what is calculated below.

		int oldValue = this->value(actor);
		int newValue = oldValue + difference;
		int totalChange = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();
			int alterValue = this->value(j);
			int change =
				std::abs(oldValue - alterValue) - std::abs(newValue - alterValue);

			if (this->lalterPopularity)
			{
				change *= pNetwork->inDegree(j);
			}

			totalChange += change;
		}

		contribution = ((double) totalChange) / this->range();

		if (this->laverage)
		{
			contribution /= pNetwork->outDegree(actor);
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
double SimilarityEffect::egoStatistic(int ego,
	double * currentValues)
{
	const Network * pNetwork = this->pNetwork();

	double statistic = 0;
	int neighborCount = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
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

	if (this->legoPopularity)
	{
		statistic *= pNetwork->inDegree(ego);
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double SimilarityEffect::egoEndowmentStatistic(int ego, const int * difference,
	double * currentValues)
{
	if (this->lalterPopularity)
	{
		throw runtime_error(string("endowmentStatistic not implemented for") +
			"average similarity x popularity alter effect and " +
			"total similarity x popularity alter effect.");
	}

	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	double similarityMean =  this->similarityMean();

	if (!this->missing(this->period(), ego) &&
		!this->missing(this->period() + 1, ego))
	{
		if (difference[ego] > 0)
		{
			if (pNetwork->outDegree(ego))
			{
				double thisStatistic = 0;

				for (IncidentTieIterator iter = pNetwork->outTies(ego);
					 iter.valid();
					 iter.next())
				{
					if (!this->missing(this->period(), iter.actor()) &&
						!this->missing(this->period() + 1, iter.actor()))
					{
						double alterValue = currentValues[iter.actor()];
						double range = this->range();
						thisStatistic += iter.value() *
							(1.0 - fabs(alterValue - currentValues[ego]) /
								range);
						thisStatistic -= similarityMean;
					}
				}

				if (this->laverage)
				{
					thisStatistic /= pNetwork->outDegree(ego);
				}

				if (this->legoPopularity)
				{
					thisStatistic *= pNetwork->inDegree(ego);
				}
				statistic = thisStatistic;

				// do the same using the current state plus difference
				// in i's value rather than current state and subtract it.
				// not sure whether this is correct.

				thisStatistic = 0;

				for (IncidentTieIterator iter = pNetwork->outTies(ego);
					 iter.valid();
					 iter.next())
				{
					if (!this->missing(this->period(), iter.actor()) &&
						!this->missing(this->period() + 1, iter.actor()))
					{
						double alterValue = currentValues[iter.actor()] +
							difference[iter.actor()];
						double range = this->range();
						thisStatistic += iter.value() *
							(1.0 - fabs(alterValue - (difference[ego] +
									currentValues[ego]))
								/ range);
						thisStatistic -= similarityMean;
					}
				}

				if (this->laverage)
				{
					thisStatistic /= pNetwork->outDegree(ego);
				}
				if (this->legoPopularity)
				{
					thisStatistic *= pNetwork->inDegree(ego);
				}

				statistic -= thisStatistic;
			}
		}
	}

	return statistic;
}

}
