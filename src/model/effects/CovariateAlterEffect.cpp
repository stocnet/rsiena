/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: https://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateAlterEffect class.
 *****************************************************************************/

#include "CovariateAlterEffect.h"
#include "network/Network.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
CovariateAlterEffect::CovariateAlterEffect(const EffectInfo * pEffectInfo,
		const bool leftThresholded, const bool rightThresholded,
		const bool squared) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lleftThresholded = leftThresholded;
	this->lrightThresholded = rightThresholded;
	this->lthreshold = pEffectInfo->internalEffectParameter();
	// to make sure that there will be no numerical equality difficulties:
	if (this->lleftThresholded)
	{
		this->lthreshold += 1e-12;
	}
	if (this->lrightThresholded)
	{
		this->lthreshold -= 1e-12;
	}
	this->lsquared = squared;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateAlterEffect::calculateContribution(int alter) const
{
	double change = 0;

	if (this->lleftThresholded)
	{
		if (this->value(alter) <= this->lthreshold)
		{
			change = 1;
		}
	}
	else
	{
		if (this->lrightThresholded)
		{
			if (this->value(alter) >= this->lthreshold)
			{
				change = 1;
			}
		}
		else
		{
			change = this->value(alter);
			if (this->lsquared)
			{
				change *= change;
			}
		}
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateAlterEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(alter))
	{
		statistic = this->calculateContribution(alter);
	}

	return statistic;
}

}
