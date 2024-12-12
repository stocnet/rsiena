/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreeContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IndegreeEffect class.
 *****************************************************************************/

#include "OutIndegreeBalanceContinuousEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/ContinuousVariable.h"
namespace siena
{

/**
 * Constructor.
 */
OutIndegreeBalanceContinuousEffect::OutIndegreeBalanceContinuousEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentContinuousEffect(pEffectInfo)
{
}


/**
 * Returns the outdegree-indegree balance of a certain actor, and thus how much this effect
 * contributes to the change in the continuous behavior.
 */
double OutIndegreeBalanceContinuousEffect::calculateChangeContribution(int ego)
{
	return this->outIndegreeBalance(ego);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double OutIndegreeBalanceContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego] * this->outIndegreeBalance(ego);
}

/**
 * Returns the outdegree-indegree balance ( O-I / O+I ) of the given actor in
 * the current network (0, if the actor has no outgoing or incoming ties).
 */
double OutIndegreeBalanceContinuousEffect::outIndegreeBalance(int ego) const
{
    double balance = this->pNetwork()->outDegree(ego) + this->pNetwork()->inDegree(ego);

    if (balance > 0)
    {
        balance = (this->pNetwork()->outDegree(ego) - this->pNetwork()->inDegree(ego)) / balance;
    }

    return balance;
}

}
