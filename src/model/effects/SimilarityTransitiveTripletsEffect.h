/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimilarityTransitiveTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * SimilarityTransitiveTripletsEffect.
 *****************************************************************************/

#ifndef SIMILARITYTRANSITIVETRIPLETSEFFECT_H_
#define SIMILARITYTRANSITIVETRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"
#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * This class defines the transitive triplets effect for non-symmetric
 * networks.
 */
class SimilarityTransitiveTripletsEffect : public CovariateDependentNetworkEffect
{
public:
	SimilarityTransitiveTripletsEffect(const EffectInfo * pEffectInfo, bool reciprocal);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the reciprocal version of the similarity transitive 
	// triplets effect is required; currently not implemented.

	bool lreciprocal {};
};

}

#endif /*SIMILARITYTRANSITIVETRIPLETSEFFECT_H_*/
