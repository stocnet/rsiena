/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CatCovariateActivityEffect.h
 *
 * Description: This file contains the definition of the
 * CatCovariateActivityEffect class.
 *****************************************************************************/

#ifndef CATCOVARIATEACTIVITYEFFECT_H_
#define CATCOVARIATEACTIVITYEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

/**
 * Same and different covariate activity effects (see manual).
 */
class CatCovariateActivityEffect : public CovariateDependentNetworkEffect
{
public:
	CatCovariateActivityEffect(const EffectInfo * pEffectInfo);
	virtual ~CatCovariateActivityEffect();
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// An array of counts of tie values
	// between ego and alters with specific covariate values
	int * lpNumberTieValues {};
	int lcovMax {};
	bool lroot {};
	double changeStat(double d, bool nosqrt) const;
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /*CATCOVARIATEACTIVITYEFFECT_H_*/
