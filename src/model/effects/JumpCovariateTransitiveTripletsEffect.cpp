/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: JumpCovariateTransitiveTripletsEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * JumpCovariateTransitiveTripletsEffect.
 *****************************************************************************/
#include <cmath>

#include "JumpCovariateTransitiveTripletsEffect.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */

JumpCovariateTransitiveTripletsEffect::JumpCovariateTransitiveTripletsEffect(
									const EffectInfo * pEffectInfo, bool reciprocal) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
// currently not used:
	{
		this->lreciprocal = reciprocal;
	}
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double JumpCovariateTransitiveTripletsEffect::calculateContribution(
                                              int alter) const
{
	// If we are introducing a tie from the ego i to the alter j, then each
	// two-path from i to j via h with v_i = v_h != v_j contributes one unit;
   // in addition, each in-star i -> h <- j with v_i = v_j != v_h
   // also contributes one unit.
   // This number is not stored in a table and is calculated from scratch.

	int contribution1 = 0;
	const Network * pNetwork = this->pNetwork();

	// The following probably can be done more efficiently 
	// using CommonNeighborIterator.
	// Iterate over ego's outgoing ties
	if (fabs(this->value(alter) - this->value(this->ego())) < EPSILON)
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
			{
				// Get the receiver of the outgoing tie.
				int h = iter.actor();
				// in 2-stars:
				if (fabs(this->value(h) - this->value(this->ego())) > EPSILON &&
				pNetwork->tieValue(alter, h) >= 1)
				{
					contribution1++ ;
				}
			}
	}
	else
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
			{
				// Get the receiver of the outgoing tie.
				int h = iter.actor();
				// 2-paths:
				if (fabs(this->value(h) - this->value(this->ego())) < EPSILON &&
				pNetwork->tieValue(h, alter) >= 1)
				{
					contribution1++ ;
				}
			}
	}
	return contribution1;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double JumpCovariateTransitiveTripletsEffect::tieStatistic(int alter)
{
	
	int statistic = 0;
	const Network * pNetwork = this->pNetwork();

	if (!this->missing(this->ego()) && !this->missing(alter) &&
		fabs(this->value(alter) - this->value(this->ego())) > EPSILON)
	{
	// The following probably can be done more efficiently 
	// using CommonNeighborIterator.
	// Iterate over ego's outgoing ties
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
			{
				// Get the receiver of the outgoing tie.
				int h = iter.actor();
				// 2-paths:
				if (fabs(this->value(h) - this->value(this->ego())) < EPSILON &&
				pNetwork->tieValue(h, alter) >= 1)
				{
					statistic++ ;
				}
			}
	}
	return statistic;
}

}
