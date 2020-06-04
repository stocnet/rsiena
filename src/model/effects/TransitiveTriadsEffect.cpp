/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveTriadsEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * TransitiveTriadsEffect.
 *****************************************************************************/

#include "TransitiveTriadsEffect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
TransitiveTriadsEffect::TransitiveTriadsEffect(
		const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveTriadsEffect::calculateContribution(int alter) const
{
	return this->pTwoPathTable()->get(alter);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double TransitiveTriadsEffect::tieStatistic(int alter)
{
	return this->pTwoPathTable()->get(alter) / 6.0;
}

}
