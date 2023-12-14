/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateAvAltEffect.h
 *
 * Description: This file contains the definition of the
 * DyadicCovariateAvAltEffect class.
 *****************************************************************************/

#ifndef DYADICCOVARIATEAVALTEFFECT_H_
#define DYADICCOVARIATEAVALTEFFECT_H_

#include "DyadicCovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 *  average alter effect weighted by dyadic covariate
 */
class DyadicCovariateAvAltEffect : public DyadicCovariateAndNetworkBehaviorEffect
{
public:
	DyadicCovariateAvAltEffect(const EffectInfo * pEffectInfo, bool divide, bool asWeight, bool outgoing);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);

protected:

private:
	// divide indicates whether there will be division by the outdegree
	bool ldivide {};
	// asWeight indicates that the dyadic covariate is used as a weight;
	// if not, used as the variable.
	bool lasWeight {};
	// lpar2 specifies that the internal effect parameter is 2
	bool lpar2 {};
	// indicates if out-ties or in-ties are to be used
	bool loutgoing {};
};

}

#endif /*DYADICCOVARIATEAVALTEFFECT_H_*/
