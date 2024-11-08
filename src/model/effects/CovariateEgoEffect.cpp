/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateEgoEffect class.
 *****************************************************************************/

#include "CovariateEgoEffect.h"
#include "model/EffectInfo.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 *
 * @param[in] pEffectInfo the effect descriptor
 */
CovariateEgoEffect::CovariateEgoEffect(const EffectInfo * pEffectInfo,
		const bool leftThresholded, const bool rightThresholded) :
	CovariateDependentNetworkEffect(pEffectInfo) {
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
}

/**
 * Constructor.
 *
 * @param[in] pEffectInfo the effect descriptor
 * @param simulatedState If `true` the value() function uses the simulated
 *        state, if any or the value at the end of the period.
 */
CovariateEgoEffect::CovariateEgoEffect(const EffectInfo * pEffectInfo,
		const bool leftThresholded, const bool rightThresholded,
		const bool simulatedState) :
	CovariateDependentNetworkEffect(pEffectInfo, simulatedState) {
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
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateEgoEffect::calculateContribution(int alter) const
{
	double contribution = 0;
	if (this->lleftThresholded)
	{
		if (this->value(this->ego()) <= this->lthreshold)
		{
			contribution = 1;
		}
	}
	else
	{
		if (this->lrightThresholded)
		{
			if (this->value(this->ego()) >= this->lthreshold)
			{
				contribution = 1;
			}
		}
		else
		{
			contribution = this->value(this->ego());
		}
	}
	return contribution;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateEgoEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(this->ego()))
	{
		statistic = this->calculateContribution(alter);
	}

	return statistic;
}


/**
 * Returns if this effect is an ego effect.
 */
bool CovariateEgoEffect::egoEffect() const
{
	return true;
}

}
