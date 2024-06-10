/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: https://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CrossCovariateActivityEffect.h
 *
 * Description: This file contains the definition of the
 * CrossCovariateActivityEffect class.
 *****************************************************************************/

#ifndef CROSSCOVARIATEACTIVITYEFFECT_H_
#define CROSSCOVARIATEACTIVITYEFFECT_H_

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
class CrossCovariateActivityEffect : public CovariateDependentNetworkEffect
{
public:
	CrossCovariateActivityEffect(const EffectInfo * pEffectInfo, bool recip);

	virtual double calculateContribution(int alter) const;	
	virtual double endowmentStatistic(Network * pLostTieNetwork);
	virtual bool egoEffect() const;

protected:
	virtual double egoStatistic(int ego, const Network * pNetwork);

private:
	bool lrecip {};
	bool lsqrt {};
	bool lthree {};
	bool lcondition1(int theAlter, double theOwnValue) const;
	double changeStat(double d, bool diffSqrt) const;
	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*CROSSCOVARIATEACTIVITYEFFECT_H_*/
