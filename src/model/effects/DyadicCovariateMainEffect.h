/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateMainEffect.h
 *
 * Description: This file contains the definition of the
 * DyadicCovariateMainEffect class.
 *****************************************************************************/

#ifndef DYADICCOVARIATEMAINEFFECT_H_
#define DYADICCOVARIATEMAINEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Dyadic covariate main effect (see manual).
 */
class DyadicCovariateMainEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	DyadicCovariateMainEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*DYADICCOVARIATEMAINEFFECT_H_*/
