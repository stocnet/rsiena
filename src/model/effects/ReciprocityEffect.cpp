/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocityEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * ReciprocityEffect.
 *****************************************************************************/

#include "ReciprocityEffect.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
ReciprocityEffect::ReciprocityEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double ReciprocityEffect::calculateContribution(int alter) const
{
	double change = 0;

	// This tie flip has an effect only if there is a tie from the alter
	// to the ego.
	if (this->inTieExists(alter))
	{
		change = 1;
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double ReciprocityEffect::tieStatistic(int alter)
{
	int statistic = 0;

	if (this->inTieExists(alter))
	{
		statistic = 1;
	}

	return statistic;
}

}
