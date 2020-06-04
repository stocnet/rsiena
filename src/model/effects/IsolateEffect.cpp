/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IsolateEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IsolateEffect class.
 *****************************************************************************/

#include <cmath>
#include "IsolateEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
IsolateEffect::IsolateEffect(const EffectInfo * pEffectInfo, bool in):
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->lin = in;
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double IsolateEffect::calculateChangeContribution(int actor,
	int difference)
{
	double value = 0;

	if (lin)
	{
	if (this->pNetwork()->inDegree(actor) == 0)
	{
		value = difference;
	}
	}
	else
	{
		if (this->pNetwork()->outDegree(actor) == 0)
		{
			value = difference;
		}
	}

	return value;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double IsolateEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;

	if (lin)
	{
	if (this->pNetwork()->inDegree(ego) == 0)
	{
		statistic = currentValues[ego];
		}
	}
	else
	{
		if (this->pNetwork()->outDegree(ego) == 0)
		{
			statistic = currentValues[ego];
		}
	}

	return statistic;
}

}
