/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveReciprocatedTripletsEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * TransitiveReciprocatedTripletsEffect.
 *****************************************************************************/

#include "TransitiveReciprocatedTripletsEffect.h"
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
TransitiveReciprocatedTripletsEffect::TransitiveReciprocatedTripletsEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveReciprocatedTripletsEffect::calculateContribution(int alter) const
{
	// If we are introducing a tie from the ego i to the alter j, then each
	// two-path from i to j with an existing tie from j to i contributes one
	// unit; in addition, each configuration i <-> h <- j also contributes one
	// unit; the number of such configurations is stored in pRBTable.

	double contribution1 = 0;

	if (this->inTieExists(alter))
	{
		contribution1 = this->pTwoPathTable()->get(alter);
	}

	return contribution1 + this->pRBTable()->get(alter);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before if there is a tie from alter to ego.
 */
double TransitiveReciprocatedTripletsEffect::tieStatistic(int alter)
{

	int statistic = 0;

	if (this->inTieExists(alter))
	{
		statistic = this->pTwoPathTable()->get(alter);
	}

	return statistic;
}

}
