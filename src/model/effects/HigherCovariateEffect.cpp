/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HigherCovariateEffect.cpp
 *
 * Description: This file contains the implementation of the
 * HigherCovariateEffect class.
 *****************************************************************************/

#include "HigherCovariateEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
HigherCovariateEffect::HigherCovariateEffect(const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double HigherCovariateEffect::calculateContribution(int alter) const
{
	double change = 0;
	double egoValue = this->value(this->ego());
	double alterValue = this->value(alter);

	if (egoValue > alterValue)
	{
		change = 1;
	}
	else if (egoValue == alterValue)
	{
		change = 0.5;
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double HigherCovariateEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (this->missing(this->ego()) || this->missing(alter))
	{
		statistic = 0.5;
	}
	else if (this->value(this->ego()) > this->value(alter))
	{
		statistic = 1;
	}
	else if (this->value(this->ego()) == this->value(alter))
	{
		statistic = 0.5;
	}

	return statistic;
}

}
