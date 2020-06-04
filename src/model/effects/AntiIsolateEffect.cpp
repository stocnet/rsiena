/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AntiIsolateEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AntiIsolateEffect class.
 *****************************************************************************/

#include "AntiIsolateEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AntiIsolateEffect::AntiIsolateEffect(
	const EffectInfo * pEffectInfo, bool outAlso, int minDegree) :
				NetworkEffect(pEffectInfo)
{
	this->lminDegree = minDegree;
	// Indicates the minimum degree threshold
	// hard-coded in EffectFactory.cpp to be 1 or 2 or 3.
	this->loutAlso = outAlso;
	// Indicates that also out-tie-isolation of alter is considered
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double AntiIsolateEffect::calculateContribution(int alter) const
{
	double change = 0;
	int degree = this->pNetwork()->inDegree(alter);

// The following could be combined in one statement
// but this would require more comparisons.
// Note that lminDegree can only be 1 or 2 or 3.
	if (lminDegree <= 1)
	{
		if ((degree <= 0) || ((degree <= 1) && (this->outTieExists(alter))))
		{
			if (loutAlso)
			{
				if (this->pNetwork()->outDegree(alter) <= 0)
				{
					change = 1;
				}
			}
			else
			{
				change = 1;
			}
		}
	}
	else // lminDegree >= 2; implies !loutAlso
	{
		if ((((degree+1) == lminDegree)  && (!this->outTieExists(alter))) ||
				((degree == lminDegree) && (this->outTieExists(alter))))
		{
			change = 1;
		}
	}

	return change;
}

/**
 * The contribution of ego to the statistic.
 * In this case, it would be better to say "contribution of alter"
 * but it is equivalent.
 * It is assumed that preprocessEgo(ego) has been called before.
 */
double AntiIsolateEffect::egoStatistic(int ego, const Network * pNetwork)
{
	double statistic = 0;

	// This is a sum over the alters; since it also applies to two-mode networks,
	// it cannot be represented as a sum over ego.
	if (ego <= 0)
	{
		for (int i = 0; i < this->pNetwork()->m(); i++)
		{
			if (this->pNetwork()->inDegree(i) >= lminDegree)
			{
				if (loutAlso) // implies a one-mode network
				{
					if (this->pNetwork()->outDegree(i) <= 0)
					{
					statistic++;
					}
				}
				else
				{
					statistic++;
				}
			}
		}
	}

	return statistic;
}


}
