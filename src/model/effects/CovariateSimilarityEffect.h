/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateSimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateSimilarityEffect class.
 *****************************************************************************/

#ifndef COVARIATESIMILARITYEFFECT_H_
#define COVARIATESIMILARITYEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-related similarity and covariate-related similarity x reciprocity
 * effects (see manual).
 */
class CovariateSimilarityEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateSimilarityEffect(const EffectInfo * pEffectInfo,
		bool reciprocal);
	CovariateSimilarityEffect(const EffectInfo * pEffectInfo,
		bool reciprocal, const bool simulatedState);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the reciprocal version of the similarity effect is
	// required

	bool lreciprocal {};
};

}

#endif /*COVARIATESIMILARITYEFFECT_H_*/
