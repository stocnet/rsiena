/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageSimilarityInDist2Effect.h
 *
 * Description: This file contains the definition of the
 * AverageSimilarityInDist2Effect class.
 *****************************************************************************/

#ifndef AVERAGESIMILARITYINDIST2EFFECT_H_
#define AVERAGESIMILARITYINDIST2EFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Total or average in-similarity distance 2 effect defined as the similarity 
 * of the ego with the average of the alter's in-neighbors.
 */
class AverageSimilarityInDist2Effect : public NetworkDependentBehaviorEffect
{
public:
	AverageSimilarityInDist2Effect(const EffectInfo * pEffectInfo,
						bool divide);
	virtual ~AverageSimilarityInDist2Effect();
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	bool ldivide {};
	// Indicates whether there will be division by the outdegree of ego
	double changesim(double zalt, double zego) const;
};

}

#endif /*AVERAGESIMILARITYINDIST2EFFECT_H_*/
