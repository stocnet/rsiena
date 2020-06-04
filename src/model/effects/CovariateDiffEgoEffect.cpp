/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDiffEgoEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDiffEgoEffect class.
 *****************************************************************************/

#include "CovariateDiffEgoEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
CovariateDiffEgoEffect::CovariateDiffEgoEffect(const EffectInfo * pEffectInfo) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateDiffEgoEffect::calculateContribution(int alter) const
{
	double change = this->value(this->ego());
	change = (this->value(alter) - change) * change;
	return change;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateDiffEgoEffect::tieStatistic(int alter)
{
	double statistic = 0;
	
	if (!(this->missing(alter) || this->missing(this->ego())))
	{
		double statistic = this->value(this->ego());
		statistic = (this->value(alter) - statistic) * statistic;
		return statistic;
	}
	
	return statistic;
}

}
