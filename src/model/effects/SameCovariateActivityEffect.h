/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: https://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateActivityEffect.h
 *
 * Description: This file contains the definition of the
 * SameCovariateActivityEffect class.
 *****************************************************************************/

#ifndef SAMECOVARIATEACTIVITYEFFECT_H_
#define SAMECOVARIATEACTIVITYEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------


/**
 * Same and different covariate activity effects (see manual).
 */
class SameCovariateActivityEffect : public CovariateDependentNetworkEffect
{
public:
	SameCovariateActivityEffect(const EffectInfo * pEffectInfo,
									bool same, bool recip);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool lsame {};
	bool lrecip {};
	bool lsqrt {};
	bool lcondition1(int theAlter, double theOwnValue) const;
	bool lcondition2(int theAlter, double theOwnValue) const;
	double changeStat(double d) const;
	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*SAMECOVARIATEACTIVITYEFFECT_H_*/
