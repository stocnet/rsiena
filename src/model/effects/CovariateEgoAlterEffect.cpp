/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateEgoAlterEffect class.
 *****************************************************************************/

#include "CovariateEgoAlterEffect.h"
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
 * @param[in] nc indicates if the non-centered version of the effect is requested
 * considered
 */
CovariateEgoAlterEffect::CovariateEgoAlterEffect(
	const EffectInfo * pEffectInfo,
	bool reciprocal, bool nc) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lreciprocal = reciprocal;
	this->lnc = nc;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateEgoAlterEffect::calculateContribution(int alter) const
{
	double change = 0;

	if (!this->lreciprocal || this->inTieExists(alter))
	{
		double egoValue = this->resolvedValue(this->ego(), this->lnc);
		double alterValue = this->resolvedValue(alter, this->lnc);
		change = egoValue * alterValue;
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateEgoAlterEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(this->ego()) && !this->missing(alter) &&
		(!this->lreciprocal || this->inTieExists(alter)))
	{
		double egoValue = this->resolvedValue(this->ego(), this->lnc);
		double alterValue = this->resolvedValue(alter, this->lnc);
		statistic = egoValue * alterValue;
	}

	return statistic;
}

}
