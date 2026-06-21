/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ActivityAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ActivityAlterEffect class.
 *****************************************************************************/

#include "ActivityAlterEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

namespace siena
{

/**
 * Constructor.
 */
ActivityAlterEffect::ActivityAlterEffect(
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
double ActivityAlterEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->averageOutDegree(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double ActivityAlterEffect::egoStatistic(int ego, double * currentValues)
{
	double egoValue = this->lnc ? currentValues[ego] + this->overallCenterMean() :
								 currentValues[ego];
	return egoValue * this->averageOutDegree(ego);
}


/**
 * Returns the total or average out-degree of the neighbors of the given actor
 * in the current network (0, if the actor has no outgoing ties).
 */
double ActivityAlterEffect::averageOutDegree(int i) const
{
	const Network * pNetwork = this->pNetwork();
	double outDegree = 0;

	if (pNetwork->outDegree(i) > 0)
	{
		for (IncidentTieIterator iter = pNetwork->outTies(i);
			iter.valid();
			iter.next())
		{
			outDegree += pNetwork->outDegree(iter.actor());
		}
		if (this->ldivide)
		{
			outDegree /= pNetwork->outDegree(i);
		}
	}

	return outDegree;
}

}
