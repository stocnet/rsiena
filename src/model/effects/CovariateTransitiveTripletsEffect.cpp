/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateTransitiveTripletsEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateTransitiveTripletsEffect.
 *****************************************************************************/

#include "CovariateTransitiveTripletsEffect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"

namespace siena
{

/**
 * Constructor.
 */
CovariateTransitiveTripletsEffect::CovariateTransitiveTripletsEffect(
	const EffectInfo * pEffectInfo) : CovariateDependentNetworkEffect(pEffectInfo)
{
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateTransitiveTripletsEffect::calculateContribution(int alter) const
{
	// If we are introducing a tie from the ego i to the alter j, then each
	// two-instar for (i,j) contributes the covariate for j; in
	// addition, each two-path i -> h -> j also contributes the covariate value for h.
	// This number is not stored in a table and is calculated from scratch.

	const Network * pNetwork = this->pNetwork();
	double contribution1 = this->value(alter) * this->pInStarTable()->get(alter);

	// The following probably can be done more efficiently using
	// CommonNeighborIterator.
	// Iterate over ego's outgoing ties
	for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
		iter.valid();
		iter.next())
	{
		// Get the receiver of the outgoing tie.
		int h = iter.actor();
		if (pNetwork->tieValue(h, alter) >= 1)
		{
			contribution1 = contribution1 + this->value(h) ;
		}
	}

	return contribution1;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateTransitiveTripletsEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(alter))
	{
		statistic = this->value(alter)
			* this->pInStarTable()->get(alter);
	}

	return statistic;
}

}
