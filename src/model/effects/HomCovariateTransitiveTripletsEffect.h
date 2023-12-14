/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HomCovariateTransitiveTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * HomCovariateTransitiveTripletsEffect.
 * Contributed by Robert Hellpap.
 *****************************************************************************/

#ifndef HOMCOVARIATETRANSITIVETRIPLETSEFFECT_H_
#define HOMCOVARIATETRANSITIVETRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"
#include "CovariateDependentNetworkEffect.h"

namespace siena
{

class HomCovariateTransitiveTripletsEffect :
                   public CovariateDependentNetworkEffect
{
public:
	HomCovariateTransitiveTripletsEffect(
                   const EffectInfo * pEffectInfo, bool reciprocal);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if a reciprocal version of the effect is required;
   // currently not implemented.

	bool lreciprocal {};
};

}

#endif /*HOMCOVARIATETRANSITIVETRIPLETSEFFECT_H_*/
