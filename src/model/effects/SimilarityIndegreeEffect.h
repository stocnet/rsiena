/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimilarityIndegreeEffect.h
 *
 * Description: This file contains the definition of the
 * SimilarityIndegreeEffect class.
 *****************************************************************************/

#ifndef SIMILARITYINDEGREEEFFECT_H_
#define SIMILARITYINDEGREEEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * This class implements several behavior effects related to indegree-similarity
 * (see manual):
 * - Average indegree-similarity
 * - Average indegree-similarity x popularity alter
 * - Total indegree-similarity
 * - Total indegree-similarity x popularity alter
 */
class SimilarityIndegreeEffect : public NetworkDependentBehaviorEffect
{
public:
	SimilarityIndegreeEffect(const EffectInfo * pEffectInfo,
		bool average,
		bool alterPopularity);
	SimilarityIndegreeEffect(const EffectInfo * pEffectInfo,
		bool average,
		bool alterPopularity, const bool simulatedState);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	bool laverage {};
	bool lalterPopularity {};
};

}

#endif /*SIMILARITYINDEGREEEFFECT_H_*/
