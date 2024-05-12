/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AllSimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * AllSimilarityEffect class.
 *****************************************************************************/

#ifndef ALLSIMILARITYEFFECT_H_
#define ALLSIMILARITYEFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

/**
 * This class implements several behavior effects related to all-similarity
 * (see manual), which do not take the network into account.
 * Constructed based on the SimilarityEffect class.
 */
class AllSimilarityEffect : public BehaviorEffect
{
public:
	AllSimilarityEffect(const EffectInfo * pEffectInfo,
		bool nearby);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	bool lnear{};
	int lp{};
};

}

#endif /*ALLSIMILARITYEFFECT_H_*/
