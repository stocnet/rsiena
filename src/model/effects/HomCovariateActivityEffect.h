/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HomCovariateActivityEffect.h
 *
 * Description: This file contains the definition of the
 * HomCovariateActivityEffect class.
 *****************************************************************************/

#ifndef HOMCOVARIATEACTIVITYEFFECT_H_
#define HOMCOVARIATEACTIVITYEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Same and different covariate activity effects (see manual).
 */
class HomCovariateActivityEffect : public CovariateDependentNetworkEffect
{
public:
	HomCovariateActivityEffect(const EffectInfo * pEffectInfo, bool same);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool lsame {};
};

}

#endif /*HOMCOVARIATEACTIVITYEFFECT_H_*/
