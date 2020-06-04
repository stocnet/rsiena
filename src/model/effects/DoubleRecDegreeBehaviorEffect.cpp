/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleRecDegreeBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DoubleRecDegreeBehaviorEffect class.
 *****************************************************************************/

#include <string>
#include <stdexcept>
#include "network/OneModeNetwork.h"
#include "DoubleRecDegreeBehaviorEffect.h"
#include "network/CommonNeighborIterator.h"
//#include "network/IncidentTieIterator.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
DoubleRecDegreeBehaviorEffect::DoubleRecDegreeBehaviorEffect(
	const EffectInfo * pEffectInfo, int secondDirection) :
		TwoNetworkDependentBehaviorEffect(pEffectInfo)
{
	if (!((secondDirection==0)||(secondDirection==1)||
									(secondDirection==2)))
	{
		throw runtime_error
			("DoubleRecDegreeBehaviorEffect: secondDirection must be 0, 1, or 2");
	}
	this->lsecondDirection = secondDirection;

	const Network * pFFirstNetwork = this->pFirstNetwork();
	const OneModeNetwork * pOFirstNetwork =
			dynamic_cast<const OneModeNetwork *>(pFFirstNetwork);
	if (!pOFirstNetwork)
	{
		throw runtime_error(
		"One-mode first network expected in DoubleRecDegreeBehaviorEffect");
	}
}

/**
 * Calculates the double degree, for use in both calculateChangeContribution
 * and egoStatistic
 */
	int DoubleRecDegreeBehaviorEffect::calculateDoubleRecDegree(int actor) const
{
	int statistic = 0;
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();

	CommonNeighborIterator iter(pFirstNetwork->outTies(actor),
								pFirstNetwork->inTies(actor));
// crashes!?

	if (this->lsecondDirection <= 0) // "F"
	{
		for ( ; iter.valid(); iter.next())
		{
			if (pSecondNetwork->tieValue(actor, iter.actor()) >= 1)
			{
				statistic++;
			}
		}
	}
	else if (this->lsecondDirection <= 1) // "B"
	{
		for ( ; iter.valid(); iter.next())
		{
			if (pSecondNetwork->tieValue(iter.actor(), actor) >= 1)
			{
				statistic++;
			}
		}
	}
	else //	if (this->lsecondDirection == 2) // "R"
	{
		for ( ; iter.valid(); iter.next())
		{
			if (pSecondNetwork->tieValue(iter.actor(), actor) >= 1 &&
					pSecondNetwork->tieValue(actor, iter.actor()) >= 1)
			{
				statistic++;
			}
		}
	}
	return statistic;
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double DoubleRecDegreeBehaviorEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->calculateDoubleRecDegree(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double DoubleRecDegreeBehaviorEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego] * this->calculateDoubleRecDegree(ego);
}

}
