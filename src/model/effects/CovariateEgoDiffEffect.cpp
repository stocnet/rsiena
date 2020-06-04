/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoDiffEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateEgoDiffEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <R_ext/Error.h>
#include "CovariateEgoDiffEffect.h"
#include "model/EffectInfo.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include <math.h>       /* floor */

namespace siena
{

/**
 * Constructor.
 *
 * @param[in] pEffectInfo the effect descriptor
 */
CovariateEgoDiffEffect::CovariateEgoDiffEffect(const EffectInfo * pEffectInfo,
					const bool plus, const bool minus):
	CovariateDependentNetworkEffect(pEffectInfo) {
	this->lplus = plus;
	this->lminus = minus;
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateEgoDiffEffect::calculateContribution(int alter) const
{
	double contribution = 0;
	int xi = floor(this->value(this->ego()));
	int deg = this->pNetwork()->outDegree(this->ego());
	if (this->outTieExists(alter))
	{
		deg--;
	}

	if (this->lplus)
	{
		if (deg >= xi)
		{
			contribution = 1;
		}
	}
	if (this->lminus)
	{
		if (deg < xi)
		{
			contribution = -1;
		}
	}

	return contribution;
}

/**
 * Calculates the statistic corresponding to the given ego. The parameter
 * pNetwork is always the current network as there are no endowment effects
 * of this kind.
 */
double CovariateEgoDiffEffect::egoStatistic(int ego,
	const Network * pNetwork)
{
	double statistic = 0;

	if (!this->missing(this->ego()))
	{
		int diffx = pNetwork->outDegree(ego) - floor(this->value(this->ego()));

		if (this->lplus)
		{
			if (diffx > 0)
			{
				statistic = diffx;
			}
		}
		if (this->lminus)
		{
			if (diffx < 0)
			{
				statistic = -diffx;
			}
		}
	}

	return statistic;
}

/**
 * Returns if this effect is an ego effect.
 */
bool CovariateEgoDiffEffect::egoEffect() const
{
	return true;
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double CovariateEgoDiffEffect::endowmentStatistic(Network * pLostTieNetwork)
{
	throw std::logic_error(
		"Positive/negative outdegree covariate difference effect: Endowment effect not supported.");
}

}
