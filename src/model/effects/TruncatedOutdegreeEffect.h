/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TruncatedOutdegreeEffect.h
 *
 * Description: This file contains the definition of the
 * TruncatedOutdegreeEffect class.
 *****************************************************************************/

#ifndef TRUNCATEDOUTDEGREEEFFECT_H_
#define TRUNCATEDOUTDEGREEEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines for right the truncated outdegree activity effect
 * defined by s_i(x) = min{x_{i+}, c}. The corresponding statistic is
 * the sum of outdegrees truncated at c over all actors.
 * For !right it defines the right-truncated outdegree activity effect
 * defined by s_i(x) = max{x_{i+}, c}.
 */
class TruncatedOutdegreeEffect : public NetworkEffect
{
public:
	TruncatedOutdegreeEffect(const EffectInfo * pEffectInfo, 
									bool right, bool outIso);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);

private:
	int lc {};
	bool lright {};
	bool lOutIso {};
};

}

#endif /*TRUNCATEDOUTDEGREEEFFECT_H_*/
