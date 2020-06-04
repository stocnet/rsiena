/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InIsolateDegreeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * InIsolateDegreeEffect class.
 *****************************************************************************/

#include "InIsolateDegreeEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "data/OneModeNetworkLongitudinalData.h"

namespace siena
{

/**
 * Constructor.
 */
InIsolateDegreeEffect::InIsolateDegreeEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double InIsolateDegreeEffect::calculateContribution(int alter) const
{
	double change = 0;
	
	if (this->pNetwork()->inDegree(this->ego()) == 0)
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
double InIsolateDegreeEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (this->pNetwork()->inDegree(this->ego()) == 0)
	{
		statistic = 1;
	}
	return statistic;
}

}
