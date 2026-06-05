/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalGwInAltDist2NCEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * TotalGwInAltDist2NCEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>

#include "TotalGwInAltDist2NCEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
TotalGwInAltDist2NCEffect::TotalGwInAltDist2NCEffect(
	const EffectInfo * pEffectInfo) :
	NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->linternalEffectParameter = pEffectInfo->internalEffectParameter();
	this->lweight = -0.01 * this->linternalEffectParameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = 1 - exp(this->lweight);

	if (this->linternalEffectParameter < 0)
	{
		throw runtime_error(
			"TotalGwInAltDist2 must have nonnegative internal effect parameter");
	}
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change their behavior by the given amount.
 */
double TotalGwInAltDist2NCEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		double weightedAlterSum = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			int alter = iter.actor();
			double total = this->totalInAlterValue(alter) +
				pNetwork->inDegree(alter) * this->overallCenterMean();
			total -= this->value(actor);
			weightedAlterSum += this->gwWeight(total);
		}

		contribution = difference * weightedAlterSum;
	}

	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * current values of the behavior variable.
 */
double TotalGwInAltDist2NCEffect::egoStatistic(int ego,
	double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		iter.valid();
		iter.next())
	{
		int alter = iter.actor();
		double total = 0;

		for (IncidentTieIterator inIter = pNetwork->inTies(alter);
			inIter.valid();
			inIter.next())
		{
			int inAlter = inIter.actor();

			if (inAlter != ego)
			{
				total += currentValues[inAlter] + this->overallCenterMean();
			}
		}

		statistic += this->gwWeight(total);
	}

	statistic *= currentValues[ego] + this->overallCenterMean();

	return statistic;
}


/**
 * Returns the geometrically weighted transform of the given alter total.
 */
double TotalGwInAltDist2NCEffect::gwWeight(double total) const
{
	if (total <= 0)
	{
		return 0;
	}

	return this->lexpmweight * (1 - pow(this->lexpfactor, total));
}

}