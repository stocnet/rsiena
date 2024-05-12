/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateTransitiveTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * SameCovariateTransitiveTripletsEffect.
 *****************************************************************************/

#ifndef SAMECOVARIATETRANSITIVETRIPLETSEFFECT_H_
#define SAMECOVARIATETRANSITIVETRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"
#include "CovariateDependentNetworkEffect.h"

namespace siena
{

class SameCovariateTransitiveTripletsEffect :
	public CovariateDependentNetworkEffect
{
public:
	SameCovariateTransitiveTripletsEffect(
			const EffectInfo * pEffectInfo, bool same);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool inequalityCondition(int a) const;
	// lsame indicates if the requirement inthe condition is "same" or "different"
	bool lsame {};
};

}

#endif /*SAMECOVARIATETRANSITIVETRIPLETSEFFECT_H_*/
