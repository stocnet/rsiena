/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AllSimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AllSimilarityEffect class.
 *****************************************************************************/
#include <cstdlib>
#include <stdexcept>
#include "AllSimilarityEffect.h"
#include "model/EffectInfo.h"
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
AllSimilarityEffect::AllSimilarityEffect(
	const EffectInfo * pEffectInfo,
	bool nearby): BehaviorEffect(pEffectInfo)
{
	this->lnear = nearby;
	if (pEffectInfo->internalEffectParameter() < 0)
	{
		throw logic_error(
			"Effect parameter for AllSimilarityEffect should be nonnegative.");
	}
	if ((pEffectInfo->internalEffectParameter() == 0) && this->lnear)
	{
		throw logic_error(
			"Effect parameter for simAllFar should be at least 1.");
	}
// Else the if-statements below do not work.
	this->lp = (pEffectInfo->internalEffectParameter());
// In C++, casting double to integers is by chopping off the decimal part.
}



/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 * A more efficient coding would be obtained by using a table of
 * the cumulative distribution of the behavior
 * (calculated before looping over all actors).
 */
double AllSimilarityEffect::calculateChangeContribution(int actor,
	int difference)
{
	int contribution = 0;
	int zi = value(actor); // this is the non-centered value
	if (difference < 0) // then change statistic for zi-1 must be taken
	{
		zi--;
	}

	for (int j = 0; j < this->n(); j++)
	{
		int zj = value(j);
		if ((zj <= zi) && (j != actor))
		{
			if (((zj > zi - this->lp) && (this->lnear)) || ((zj <= zi - this->lp) && (!this->lnear)))
			{
				contribution--;
			}
		}
		if ((zj > zi) && (j != actor))
		{
			if (((zj <= zi + this->lp) && (this->lnear)) || ((zj > zi + this->lp) && (!this->lnear)))
			{
				contribution++;
			}
		}
	}
	return difference * contribution;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AllSimilarityEffect::egoStatistic(int ego,
	double * currentValues)
{
	double statistic = 0;
	double zi = currentValues[ego]; // this is the non-centered value

	for (int j = 0; j < this->n(); j++)
	{
		double dif = currentValues[j]-zi;
		if ((dif < 0) && (j != ego))
		{
			if (((dif > - this->lp) && (this->lnear)) || ((dif <= - this->lp) && (!this->lnear)))
			{
				statistic = statistic + this->lp + dif;
			}
		}
		if ((dif > 0) && (j != ego))
		{
			if (((dif < this->lp) && (this->lnear)) || ((dif > this->lp) && (!this->lnear)))
			{
				statistic = statistic + this->lp - dif;
			}
		}
		if ((dif == 0) && (j != ego) && (this->lnear))
		{
			statistic = statistic + this->lp;
		}
	}
	return statistic;
}
}
