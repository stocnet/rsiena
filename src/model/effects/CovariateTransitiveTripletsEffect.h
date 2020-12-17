/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateTransitiveTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * CovariateTransitiveTripletsEffect.
 *****************************************************************************/

#ifndef COVARIATETRANSITIVETRIPLETSEFFECT_H_
#define COVARIATETRANSITIVETRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"
#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * This class defines the covariate transitive triplets effect for non-symmetric
 * networks.
 */
class CovariateTransitiveTripletsEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateTransitiveTripletsEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

};

}

#endif /*COVARIATETRANSITIVETRIPLETSEFFECT_H_*/
