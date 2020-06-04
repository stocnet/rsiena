/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TruncatedOutdegreeEffect2.h
 *
 * Description: This file contains the definition of the
 * TruncatedOutdegreeEffect2 class.
 *****************************************************************************/

#ifndef TRUNCATEDOUTDEGREEEFFECT2_H_
#define TRUNCATEDOUTDEGREEEFFECT2_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the outdegree activity effect defined by
 * s_i(x) = min{x_{i+}, c}. The corresponding statistic is
 * the sum of outdegrees truncated at c over all actors.
 */
class TruncatedOutdegreeEffect2 : public NetworkEffect
{
public:
	TruncatedOutdegreeEffect2(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;
protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);

private:
	double lc;
};

}

#endif /*TRUNCATEDOUTDEGREEEFFECT2_H_*/
