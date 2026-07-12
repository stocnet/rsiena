/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PopularityAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * PopularityAlterEffect class.
 *****************************************************************************/

#include "PopularityAlterEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

namespace siena
{

/**
 * Constructor.
 */
PopularityAlterEffect::PopularityAlterEffect(
		const EffectInfo * pEffectInfo, bool divide, bool nc) :
	NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
	this->lnc = nc;
	// Indicates whether behavior values will not be centered
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double PopularityAlterEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->averageInDegree(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double PopularityAlterEffect::egoStatistic(int ego, double * currentValues)
{
	double egoValue = this->lnc ? currentValues[ego] + this->overallCenterMean() :
								 currentValues[ego];
	return egoValue * this->averageInDegree(ego);
}


/**
 * Returns the total or average in-degree of the neighbors of the given actor in
 * the current network (0, if the actor has no outgoing ties).
 */
double PopularityAlterEffect::averageInDegree(int i) const
{
	const Network * pNetwork = this->pNetwork();
	double inDegree = 0;

	if (pNetwork->outDegree(i) > 0)
	{
		for (IncidentTieIterator iter = pNetwork->outTies(i);
			iter.valid();
			iter.next())
		{
			inDegree += pNetwork->inDegree(iter.actor());
		}
		if (this->ldivide)
		{
		inDegree /= pNetwork->outDegree(i);
		}
	}

	return inDegree;
}

}
