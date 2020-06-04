/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateReciprocityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DyadicCovariateReciprocityEffect class.
 *****************************************************************************/

#include "DyadicCovariateReciprocityEffect.h"
#include "network/Network.h"
#include "network/CommonNeighborIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
DyadicCovariateReciprocityEffect::DyadicCovariateReciprocityEffect(
	const EffectInfo * pEffectInfo) :
		DyadicCovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double DyadicCovariateReciprocityEffect::calculateContribution(
	int alter) const
{
	double change = 0;
	int ego = this->ego();
	
	if (this->inTieExists(alter))
	{
		if (this->constantCovariate() && !this->missing(ego, alter))
		{
			change = this->value(ego, alter);
		}
		else
		{
			change = this->value(ego, alter);
		}
	}
		
	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double DyadicCovariateReciprocityEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (this->inTieExists(alter) && !this->missing(this->ego(), alter))
	{
		statistic = this->value(this->ego(), alter);
	}

	return statistic;
}

}
