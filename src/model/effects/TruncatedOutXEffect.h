/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: https://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TruncatedOutXEffect.h
 *
 * Description: This file contains the definition of the
 * TruncatedOutXEffect class.
 *****************************************************************************/

#ifndef TRUNCATEDOUTXEFFECT_H_
#define TRUNCATEDOUTXEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * This class defines the right-truncated outdegree activity effect
 * defined by s_i(x) = max{\sum_j (x_{ij} I{v_j > 0}) - p, 0}.
 */
class TruncatedOutXEffect : public CovariateDependentNetworkEffect
{
public:
	TruncatedOutXEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);

private:
	int lc {};
};

}

#endif /*TRUNCATEDOUTXEFFECT_H_*/
