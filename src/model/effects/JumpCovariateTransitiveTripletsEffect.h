/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: JumpCovariateTransitiveTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * JumpCovariateTransitiveTripletsEffect.
 *****************************************************************************/

#ifndef JUMPCOVARIATETRANSITIVETRIPLETSEFFECT_H_
#define JUMPCOVARIATETRANSITIVETRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"
#include "CovariateDependentNetworkEffect.h"

namespace siena
{

class JumpCovariateTransitiveTripletsEffect : 
                   public CovariateDependentNetworkEffect
{
public:
	JumpCovariateTransitiveTripletsEffect(
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

#endif /*JUMPCOVARIATETRANSITIVETRIPLETSEFFECT_H_*/
