/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AlterCovariateActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AlterCovariateActivityEffect class.
 *****************************************************************************/

#include "AlterCovariateActivityEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
AlterCovariateActivityEffect::AlterCovariateActivityEffect(
		const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double AlterCovariateActivityEffect::calculateContribution(int alter) const
{
	double altervalue = this->value(alter);
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();
	for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
	{
		contribution += this->value(iter.actor());
	}
	if (this->outTieExists(alter))
	{
		contribution -= altervalue;
	}
	contribution *= (2 * altervalue);
	contribution += (altervalue * altervalue);
	return contribution;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double AlterCovariateActivityEffect::tieStatistic(int alter)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (!(this->missing(alter)))
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
				iter.valid();
				iter.next())
		{
			// Get the receiver of the outgoing tie.
			int h = iter.actor();
			if (!this->missing(h))
			{
				contribution += this->value(h);
			}
		}
		contribution *= this->value(alter);
	}
	return contribution;
}

}
