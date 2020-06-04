/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AlterCovariateActivityEffect.h
 *
 * Description: This file contains the definition of the
 * AlterCovariateActivityEffect class.
 *****************************************************************************/

#ifndef ALTERCOVARIATEACTIVITYEFFECT_H_
#define ALTERCOVARIATEACTIVITYEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Alter covariate activity effects (see manual).
 */
class AlterCovariateActivityEffect : public CovariateDependentNetworkEffect
{
public:
	AlterCovariateActivityEffect(const EffectInfo * pEffectInfo);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*ALTERCOVARIATEACTIVITYEFFECT_H_*/
