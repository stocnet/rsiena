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

#include "SameCovariateTransitiveReciprocatedTripletsEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

#include <cstdlib>

namespace siena
{

/**
 * Constructor.
 */

SameCovariateTransitiveReciprocatedTripletsEffect::SameCovariateTransitiveReciprocatedTripletsEffect(
		const EffectInfo * pEffectInfo, bool same):
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsame = same;
}

bool SameCovariateTransitiveReciprocatedTripletsEffect::inequalityCondition(int a) const
{
	if (lsame)
	{
		return (abs(a) < EPSILON);
	}
	else
	{
		return (abs(a) >= EPSILON);
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double SameCovariateTransitiveReciprocatedTripletsEffect::calculateContribution(int alter) const
{
	// If we are introducing a tie from the ego i to the alter j, then each
	// two-path from i to j with an existing tie from j to i contributes one
	// unit; in addition, each configuration i <-> h <- j also contributes one
	// unit; the number of such configurations is stored in pRBTable.

	double contribution1 = 0;
	const Network * pNetwork = this->pNetwork();
	const OneModeNetwork * pOneModeNetwork =
		dynamic_cast<const OneModeNetwork *>(pNetwork);

	if (this->inTieExists(alter) &&
				this->inequalityCondition(this->value(alter) - this->value(this->ego())))
	{
		contribution1 = this->pTwoPathTable()->get(alter);
	}

	for (CommonNeighborIterator iter = pOneModeNetwork->reciprocatedTies(this->ego());
		iter.valid();
		iter.next())
	{
		int j = iter.actor();
		if ((j != alter) &&
				this->inequalityCondition(this->value(j) - this->value(this->ego()))
					&& (pNetwork->tieValue(j, alter) >= 1))
		{
			contribution1++;
		}
	}
	return contribution1;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before if there is a tie from alter to ego.
 */
double SameCovariateTransitiveReciprocatedTripletsEffect::tieStatistic(int alter)
{

	int statistic = 0;

	if (this->inTieExists(alter) && !this->missing(this->ego()) && !this->missing(alter)
			&& this->inequalityCondition(this->value(alter) - this->value(this->ego())))
	{
		statistic = this->pTwoPathTable()->get(alter);
	}

	return statistic;
}

}
