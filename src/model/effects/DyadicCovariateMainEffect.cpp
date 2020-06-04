/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateMainEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DyadicCovariateMainEffect class.
 *****************************************************************************/
#include <R_ext/Print.h>
#include "DyadicCovariateMainEffect.h"
#include "network/Network.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
DyadicCovariateMainEffect::DyadicCovariateMainEffect(
	const EffectInfo * pEffectInfo) :
		DyadicCovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double DyadicCovariateMainEffect::calculateContribution(int alter) const
{
	double change = 0;
	int ego = this->ego();

	// Exclude missings for Constant Covariates rather than replace them by
	// the mean. For changing covariates we need difference between period
	// and global means here. For constant ones we need to overwrite by 0 after
	// centering, which is equivalent to excluding.
	if (this->constantCovariate())
	{
		if (!this->missing(ego, alter))
		{
			change = this->value(ego, alter);
			//	Rprintf("const %d %d %f \n", ego, alter, change);
		}
	}
	else
	{
		change = this->value(ego, alter);
		//	Rprintf("changing %d %d %f \n", ego, alter, change);
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double DyadicCovariateMainEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(this->ego(), alter))
	{
		statistic = this->value(this->ego(), alter);
		//	Rprintf("stats %d %d %f \n", this->ego(), alter, statistic);
	}

	return statistic;
}

}
