/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveReciprocatedTriplets2Effect.cpp
 *
 * Description: This file contains the implementation of the class
 * TransitiveReciprocatedTriplets2Effect.
 *****************************************************************************/

#include "TransitiveReciprocatedTriplets2Effect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
TransitiveReciprocatedTriplets2Effect::TransitiveReciprocatedTriplets2Effect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveReciprocatedTriplets2Effect::calculateContribution(int alter) const
{
	// If we are introducing a tie from the ego i to the alter j, then each
	// configuration i -> h <-> j contributes one unit; the number of such
	// configurations is stored in pFRTable.
	return this->pFRTable()->get(alter);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before if there is a tie from alter to ego.
 */
double TransitiveReciprocatedTriplets2Effect::tieStatistic(int alter)
{
	return this->pFRTable()->get(alter);
}

}
