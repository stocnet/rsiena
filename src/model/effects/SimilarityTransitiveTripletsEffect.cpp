/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimilarityTransitiveTripletsEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * SimilarityTransitiveTripletsEffect.
 *****************************************************************************/

#include "SimilarityTransitiveTripletsEffect.h"
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
SimilarityTransitiveTripletsEffect::SimilarityTransitiveTripletsEffect(
	const EffectInfo * pEffectInfo,
	bool reciprocal) : CovariateDependentNetworkEffect(pEffectInfo)
{
	// currently not used:
	{
		this->lreciprocal = reciprocal;
	}
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double SimilarityTransitiveTripletsEffect::calculateContribution(int alter) const
{
	// If we are introducing a tie from the ego i to the alter j, then each
	// two-path from i to j with v_i = v_j contributes the similarity i-j; in
	// addition, each in-star i -> h <- j also contributes the similarity i-h.
	// This number is not stored in a table and is calculated from scratch.

	const Network * pNetwork = this->pNetwork();
	double contribution1 = this->similarity(this->ego(), alter)
		* this->pTwoPathTable()->get(alter);

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
			contribution1 = contribution1 + this->similarity(this->ego(), h) ;
		}
	}

	return contribution1;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double SimilarityTransitiveTripletsEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(this->ego()) && !this->missing(alter))
	{
		statistic = this->similarity(this->ego(), alter)
			* this->pTwoPathTable()->get(alter);
	}

	return statistic;
}

}
