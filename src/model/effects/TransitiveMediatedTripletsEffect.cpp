/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveMediatedTripletsEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * TransitiveMediatedTripletsEffect.
 *****************************************************************************/

#include "TransitiveMediatedTripletsEffect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
TransitiveMediatedTripletsEffect::TransitiveMediatedTripletsEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveMediatedTripletsEffect::calculateContribution(
	int alter) const
{
	// The number of out-stars to the ego and the given alter is the amount
	// of change for this tie flip.

	return this->pOutStarTable()->get(alter);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double TransitiveMediatedTripletsEffect::tieStatistic(int alter)
{
	return this->pOutStarTable()->get(alter);
}

}
