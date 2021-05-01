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
#include "model/EffectInfo.h"
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
	bool egoPopularity,
	bool hi,
	bool lo) : NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->laverage = average;
	this->lalterPopularity = alterPopularity;
	this->legoPopularity = egoPopularity;
	this->lhi = hi;
	this->llo = lo;
	this->lcenter = true;
}

/**
 * Constructor.
 * @param[in] average indicates if one of the average effects is required
 * @param[in] alterPopularity indicates if the similarity scores have to
 * be multiplied by the in-degrees of alters
 * @param simulatedState If `true` the value() function uses the simulated
 *        state, if any or the value at the end of the period.
 */
SimilarityEffect::SimilarityEffect(
	const EffectInfo * pEffectInfo,
	bool average,
	bool alterPopularity,
	bool egoPopularity,
	bool hi,
	bool lo,
	const bool simulatedState) :
		NetworkDependentBehaviorEffect(pEffectInfo,simulatedState)
{
	this->laverage = average;
	this->lalterPopularity = alterPopularity;
	this->legoPopularity = egoPopularity;
	this->lhi = hi;
	this->llo = lo;
	this->lcenter = (pEffectInfo->internalEffectParameter() <= 1);
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
		// s_i(x) = avg(sim(v_i, v_j) - similarityMean) over all neighbors
		// j of i.
		// sim(v_i, v_j) = 1.0 - |v_i - v_j| / observedRange

		int totalChange = 0;

		if (this->lalterPopularity)
		{
			if (difference > 0)
			{
				if (this->lhi)
				{
					totalChange	+= this->numberAlterHigherPop(actor);
				}
				if (this->llo)
				{
					totalChange -= this->numberAlterEqualPop(actor);
					totalChange -= this->numberAlterLowerPop(actor);
				}
			}
			else if (difference < 0)
			{
				if (this->lhi)
				{
					totalChange -= this->numberAlterHigherPop(actor);
					totalChange -= this->numberAlterEqualPop(actor);
				}
				if (this->llo)
				{
					totalChange += this->numberAlterLowerPop(actor);
				}
			}
		}
		else
		{
			if (difference > 0)
			{
				if (this->lhi)
				{
					totalChange	+= this->numberAlterHigher(actor);
				}
				if (this->llo)
				{
					totalChange -= this->numberAlterEqual(actor);
					totalChange -= this->numberAlterLower(actor);
				}
			}
			else if (difference < 0)
			{
				if (this->lhi)
				{
					totalChange -= this->numberAlterHigher(actor);
					totalChange -= this->numberAlterEqual(actor);
				}
				if (this->llo)
				{
					totalChange += this->numberAlterLower(actor);
				}
			}
		}

		contribution = ((double) totalChange) / this->range();

		if (this->laverage)
		{
			contribution /= pNetwork->outDegree(actor);
		}
		else
		{
			if (this->lhi && this->llo && this->lcenter)
			{
				contribution -= (pNetwork->outDegree(actor) * this->similarityMean());
			}
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
	double statistic = 0;
	int neighborCount = 0;
	int totalCount = 0;
	const Network * pNetwork = this->pNetwork();

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j))
		{
			double difValues = currentValues[j] - currentValues[ego];
			if (this->lalterPopularity)
			{
				difValues *= pNetwork->inDegree(j);
			}

			if (this->lhi)
			{
				if (difValues > 0)
				{
					statistic += difValues;
				}
			}
			if (this->llo)
			{
				if (difValues < 0)
				{
					statistic -= difValues;
				}
			}
			neighborCount++;
			if (this->lalterPopularity)
			{
				totalCount += pNetwork->inDegree(j);
			}
		}
	} // until here loop outTies j

	if (!this->lalterPopularity)
	{
		totalCount = neighborCount;
	}
			
	statistic = totalCount - (statistic / this->range());

// Center similarity for the avSim and related effects
	if (this->lhi && this->llo && this->lcenter)
	{
		statistic = statistic - (totalCount * this->similarityMean());
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
	int neighborCount = 0;
	const Network * pNetwork = this->pNetwork();

	if (!this->missing(this->period(), ego) &&
		!this->missing(this->period() + 1, ego))
	{
		if (difference[ego] > 0)
		{
			if (pNetwork->outDegree(ego))
			{
				double thisStatistic = 0;
				double egoValue = currentValues[ego];

				for (IncidentTieIterator iter = pNetwork->outTies(ego);
					 iter.valid();
					 iter.next())
				{
					if (!this->missing(this->period(), iter.actor()) &&
						!this->missing(this->period() + 1, iter.actor()))
					{
						int j = iter.actor();
						double difValues = currentValues[j] - egoValue;
						if (this->lhi)
						{
							if (difValues > 0)
							{
								thisStatistic += difValues;
							}
						}
						if (this->llo)
						{
							if (difValues < 0)
							{
								thisStatistic -= difValues;
							}
						}
						neighborCount++;
					}
				}
				statistic = thisStatistic;

				// do the same using the current state plus difference
				// in i's value rather than current state and subtract it.
				// not sure whether this is correct.

				thisStatistic = 0;

				egoValue = currentValues[ego] + difference[ego];

				for (IncidentTieIterator iter = pNetwork->outTies(ego);
					 iter.valid();
					 iter.next())
				{
					int j = iter.actor();
					if (!this->missing(this->period(), j) &&
						!this->missing(this->period() + 1, j))
					{
						double difValues = currentValues[j] + difference[j] - egoValue;

						if (this->lhi)
						{
							if (difValues > 0)
							{
								thisStatistic += difValues;
							}
						}
						if (this->llo)
						{
							if (difValues < 0)
							{
								thisStatistic -= difValues;
							}
						}
					}
				}
				statistic -= thisStatistic;
				if (this->laverage && neighborCount > 0)
				{
					statistic /= neighborCount;
				}
			}
		}
	}

	return statistic;
}

}
