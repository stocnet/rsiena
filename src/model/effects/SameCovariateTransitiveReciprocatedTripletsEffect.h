/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateTransitiveReciprocatedTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * SameCovariateTransitiveReciprocatedTripletsEffect.
 *****************************************************************************/

#ifndef SAMECOVARIATETRANSITIVERECIPROCATEDTRIPLETSEFFECT_H_
#define SAMECOVARIATETRANSITIVERECIPROCATEDTRIPLETSEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * This class defines the reciprocated transitive triplets effect for non-symmetric
 * networks.
 */
class SameCovariateTransitiveReciprocatedTripletsEffect : public CovariateDependentNetworkEffect
{
public:
	SameCovariateTransitiveReciprocatedTripletsEffect(const EffectInfo * pEffectInfo, bool same);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool inequalityCondition(int a) const;
	// lsame indicates if the requirement inthe condition is "same" or "different"
	bool lsame {};
};

}

#endif /*SAMECOVARIATETRANSITIVERECIPROCATEDTRIPLETSEFFECT_H_*/
