/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoSquaredEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateEgoSquaredEffect class.
 *****************************************************************************/

#include "CovariateEgoSquaredEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
CovariateEgoSquaredEffect::CovariateEgoSquaredEffect(const EffectInfo * pEffectInfo) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateEgoSquaredEffect::calculateContribution(int alter) const
{
	double change = this->value(this->ego());
	change *= change;
	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateEgoSquaredEffect::tieStatistic(int alter)
{
	double statistic = 0;
	if (!this->missing(this->ego()))
	{
		statistic = this->value(this->ego());
		statistic *= statistic;
	}
	return statistic;
}

/**
 * Returns if this effect is an ego effect.
 */
bool CovariateEgoSquaredEffect::egoEffect() const
{
	return true;
}

}
