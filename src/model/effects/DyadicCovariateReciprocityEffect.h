/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateReciprocityEffect.h
 *
 * Description: This file contains the definition of the
 * DyadicCovariateReciprocityEffect class.
 *****************************************************************************/

#ifndef DYADICCOVARIATERECIPROCITYEFFECT_H_
#define DYADICCOVARIATERECIPROCITYEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Dyadic covariate x reciprocity effect (see manual).
 */
class DyadicCovariateReciprocityEffect :
	public DyadicCovariateDependentNetworkEffect
{
public:
	DyadicCovariateReciprocityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*DYADICCOVARIATERECIPROCITYEFFECT_H_*/
