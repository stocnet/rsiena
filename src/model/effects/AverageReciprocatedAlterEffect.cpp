/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageReciprocatedAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageReciprocatedAlterEffect class.
 *****************************************************************************/

#include <string>
#include <stdexcept>
#include "AverageReciprocatedAlterEffect.h"
#include "network/CommonNeighborIterator.h"
#include "network/OneModeNetwork.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AverageReciprocatedAlterEffect::AverageReciprocatedAlterEffect(
	const EffectInfo * pEffectInfo, bool divide) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AverageReciprocatedAlterEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(string("One-mode network expected in ") +
			"AverageReciprocatedAlterEffect");
	}

	if (pNetwork->reciprocalDegree(actor) > 0)
	{
		double totalAlterValue = 0;

		for (CommonNeighborIterator iter = pNetwork->reciprocatedTies(actor);
			iter.valid();
			iter.next())
		{
			double alterValue = this->centeredValue(iter.actor());
			totalAlterValue += alterValue;
		}

		contribution = difference * totalAlterValue;
		if (this->ldivide)
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
double AverageReciprocatedAlterEffect::egoStatistic(int i,
	double * currentValues)
{
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(string("One-mode network expected in ") +
			"AverageReciprocatedAlterEffect");
	}

	double statistic = 0;
	int neighborCount = 0;

	for (CommonNeighborIterator iter = pNetwork->reciprocatedTies(i);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j))
		{
			statistic += currentValues[j];
			neighborCount++;
		}
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

}
