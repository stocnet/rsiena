/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleDegreeBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DoubleDegreeBehaviorEffect class.
 *****************************************************************************/

#include <string>
#include <stdexcept>
#include "DoubleDegreeBehaviorEffect.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
DoubleDegreeBehaviorEffect::DoubleDegreeBehaviorEffect(
	const EffectInfo * pEffectInfo, bool firstDirection, int secondDirection) :
		TwoNetworkDependentBehaviorEffect(pEffectInfo)
{
	if (!((secondDirection==0)||(secondDirection==1)||
									(secondDirection==2)))
	{
		throw runtime_error
			("DoubleDegreeBehaviorEffect: secondDirection must be 0, 1, or 2");
	}
	this->lfirstDirection = firstDirection;
	this->lsecondDirection = secondDirection;
	this->lsubtract = (pEffectInfo->internalEffectParameter() >= 2);

//	if (secondDirection == 2)
//	{
//		const OneModeNetwork * pSSecondNetwork =
//			dynamic_cast<const OneModeNetwork *>(this->pSecondNetwork()); this crashes
//		if (!pSSecondNetwork)
//		{
//			throw runtime_error(
//			"One-mode second network expected in DoubleDegreeBehaviorEffect");
//		}
//	}
}

/**
 * Calculates the double degree, for use in both calculateChangeContribution
 * and egoStatistic
 */
	int DoubleDegreeBehaviorEffect::calculateDoubleDegree(int actor) const
{
	int statistic = 0;
	IncidentTieIterator iter;

	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();

	if (this->lfirstDirection) // "F"
	{
		iter = pFirstNetwork->outTies(actor);
	}
	else // "B"
	{
		iter = pFirstNetwork->inTies(actor);
	}

	if (this->lsecondDirection <= 0) // "F"
	{
		for (  ; iter.valid(); iter.next())
		{
			if (pSecondNetwork->tieValue(actor, iter.actor()) >= 1)
			{
				statistic++;
			}
		}
	}
	else if (this->lsecondDirection <= 1) // "B"
	{
		for (  ; iter.valid(); iter.next())
		{
			if (pSecondNetwork->tieValue(iter.actor(), actor) >= 1)
			{
				statistic++;
			}
		}
	}
	else // (lsecondDirection == 2) // "R"
	{
		for (  ; iter.valid(); iter.next())
		{
			if ((pSecondNetwork->tieValue(iter.actor(), actor) >= 1)
				&& (pSecondNetwork->tieValue(actor, iter.actor()) >= 1))
			{
				statistic++;
			}
		}
	}
	if (this->lsubtract)
	{		
		if (this->lfirstDirection) // "F"
		{
			statistic -= pFirstNetwork->outDegree(actor);
		}
		else // "B"
		{
			statistic -= pFirstNetwork->inDegree(actor);
		}
	}
	return statistic;
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double DoubleDegreeBehaviorEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->calculateDoubleDegree(actor);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double DoubleDegreeBehaviorEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego] * this->calculateDoubleDegree(ego);
}

}
