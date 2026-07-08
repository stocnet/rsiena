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
	const EffectInfo * pEffectInfo, bool forward) :
	NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->lparameter = pEffectInfo->internalEffectParameter();
	this->lweight = -0.01 * this->lparameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = 1 - exp(this->lweight);
	this->lforward = forward;

	if (this->lparameter < 0)
	{
		throw runtime_error(
			"TotalGwInAltDist2_nc must have nonnegative internal effect parameter");
	}
}


/**
 * Initializes this effect, then builds the precomputed geometrically
 * weighted transform lookup. The "total" evaluated in
 * calculateChangeContribution()/egoStatistic() is a sum of raw (non-
 * centered) int-valued behavior scores of at most n() actors, each
 * <= max(); n() * max() is therefore the tight upper bound on the index,
 * so sizing to n() * max() + 1 makes every reachable index valid 
 * The value-sum index is only well defined for
 * non-negative behavior (a negative value could push it below zero),
 * matching the precondition already documented for the analogous
 * diffusion rate effects.
 */
void TotalGwInAltDist2NCEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkDependentBehaviorEffect::initialize(pData, pState, period, pCache);

	if (this->min() < 0)
	{
		throw runtime_error(
			"TotalGwInAltDist2_nc requires a non-negative behavior variable");
	}

	double pow_ = 1;
	int size = this->pNetwork()->n() * this->max() + 1;
	this->lcumulativeWeight.resize(size); // default values 0
	for (int i = 1; i < size; i++)
	{
		pow_ *= this->lexpfactor;
		this->lcumulativeWeight[i] = this->lexpmweight * (1 - pow_);
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
			int degree = this->lforward ? pNetwork->outDegree(alter) :
											 pNetwork->inDegree(alter);
			double altervalue = this->lforward ? this->totalAlterValue(alter) :
												 this->totalInAlterValue(alter);
			double total = altervalue + degree * this->overallCenterMean();
			if (this->lforward)
			{
				// totalAlterValue(alter) sums over alter's OUT-neighbors, so
				// actor is only a member of that sum if alter has a tie back
				// to actor; not guaranteed just because actor -> alter exists.
				if (pNetwork->tieValue(alter, actor) == 1)
				{
					total -= this->value(actor);
				}
			}
			else
			{
				// totalInAlterValue(alter) sums over alter's IN-neighbors, and
				// actor is always one of them here: alter was only reached via
				// actor's own out-tie to it, i.e. there always is a tie
				// actor -> alter.
				total -= this->value(actor);
			}
			weightedAlterSum += this->lcumulativeWeight[(int) total];
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

		for (IncidentTieIterator iter = this->lforward ?
				pNetwork->outTies(alter) : pNetwork->inTies(alter);
			iter.valid();
			iter.next())
		{
			int alter = iter.actor();

			if (alter != ego)
			{
				total += currentValues[alter] + this->overallCenterMean();
			}
		}

		statistic += this->lcumulativeWeight[(int) total];
	}

	statistic *= currentValues[ego] + this->overallCenterMean();

	return statistic;
}

}