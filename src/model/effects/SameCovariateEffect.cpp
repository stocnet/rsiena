/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateEffect.cpp
 *
 * Description: This file contains the implementation of the
 * SameCovariateEffect class.
 *****************************************************************************/

#include <cmath>

#include "SameCovariateEffect.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] reciprocal indicates if only reciprocal ties have to be
 * considered
 */
SameCovariateEffect::SameCovariateEffect(const EffectInfo * pEffectInfo,
		bool reciprocal) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lreciprocal = reciprocal;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double SameCovariateEffect::calculateContribution(int alter) const
{
	int change = 0;

	if (!this->lreciprocal || this->inTieExists(alter))
	{
		if (fabs(this->value(alter) - this->value(this->ego())) < EPSILON)
		{
			change = 1;
		}
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double SameCovariateEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(this->ego()) && !this->missing(alter) &&
		(!this->lreciprocal || this->inTieExists(alter)) &&
		fabs(this->value(alter) - this->value(this->ego())) < EPSILON)
	{
		statistic = 1;
	}

	return statistic;
}

}
