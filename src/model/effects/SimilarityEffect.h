/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * SimilarityEffect class.
 *****************************************************************************/

#ifndef SIMILARITYEFFECT_H_
#define SIMILARITYEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * This class implements several behavior effects related to similarity
 * (see manual):
 * - Average similarity
 * - Average similarity x popularity alter
 * - Total similarity
 * - Total similarity x popularity alter
 * - Average similarity x popularity ego
 */
class SimilarityEffect : public NetworkDependentBehaviorEffect
{
public:
	SimilarityEffect(const EffectInfo * pEffectInfo,
		bool average,
		bool alterPopularity,
		bool egoPopularity,
		bool hi,
		bool lo);
	SimilarityEffect(const EffectInfo * pEffectInfo,
		bool average,
		bool alterPopularity,
		bool egoPopularity,
		bool hi,
		bool lo, const bool simulatedState);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	bool laverage {};
	bool lalterPopularity {};
	bool legoPopularity {};
	bool lhi {};
	bool llo {};
	bool lcenter {};
};

}

#endif /*SIMILARITYEFFECT_H_*/
