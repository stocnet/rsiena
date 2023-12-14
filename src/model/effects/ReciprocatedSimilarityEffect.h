/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocatedSimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeEffect class.
 *****************************************************************************/

#ifndef RECIPROCATEDSIMILARITYEFFECT_H_
#define RECIPROCATEDSIMILARITYEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * This class implements several behavior effects related to similarity
 * and reciprocity (see manual):
 * - Average similarity x reciprocity
 * - Average similarity x popularity alter x reciprocity
 * - Total similarity x reciprocity
 * - Total similarity x popularity alter x reciprocity
 */
class ReciprocatedSimilarityEffect : public NetworkDependentBehaviorEffect
{
public:
	ReciprocatedSimilarityEffect(const EffectInfo * pEffectInfo,
		bool average,
		bool alterPopularity);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	bool laverage {};
	bool lalterPopularity {};
};

}

#endif /*RECIPROCATEDSIMILARITYEFFECT_H_*/

