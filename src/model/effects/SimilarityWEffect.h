/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimilarityWEffect.h
 *
 * Description: This file contains the definition of the
 * SimilarityWEffect class.
 *****************************************************************************/

#ifndef SIMILARITYWEFFECT_H_
#define SIMILARITYWEFFECT_H_

#include "DyadicCovariateAndNetworkBehaviorEffect.h"

namespace siena
{

/**
 * This class implements several behavior effects related to similarity weighted by a dyadic covariate W
 * (see manual):
 * - Average similarity weighted by W
 * - Average similarity x popularity alter weighted by W
 * - Total similarity weighted by W
 * - Total similarity x popularity alter weighted by W
 * - Average similarity x popularity ego weighted by W
 * - Currently, the interactions with popularity are not in the EffectFactory.
 * - if those are going to be implemented, note what must be done with totalWeightValue
 */
class SimilarityWEffect : public DyadicCovariateAndNetworkBehaviorEffect
{
public:
	SimilarityWEffect(const EffectInfo * pEffectInfo,
		bool average,
		bool alterPopularity,
		bool egoPopularity);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	bool laverage {};
	bool lalterPopularity {};
	bool legoPopularity {};
	// lpar2 specifies that the internal effect parameter is 2
	bool lpar2 {};
};

}

#endif /*SIMILARITYWEFFECT_H_*/
