/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAvAltEffect.h
 *
 * Description: This file contains the definition of the
 * AltersCovariateAvAltEffect class.
 *****************************************************************************/

#ifndef ALTERSCOVARIATEAVALTEFFECT_H_
#define ALTERSCOVARIATEAVALTEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Alters covariate average similarity effect (not in manual)
 */
class AltersCovariateAvAltEffect :
public CovariateAndNetworkBehaviorEffect
{
public:
	AltersCovariateAvAltEffect(const EffectInfo * pEffectInfo, bool divide);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

protected:

private:
	// divide indicates whether there will be division by the outdegree
	bool ldivide{};
};

}

#endif /*ALTERSCOVARIATEAVALTEFFECT_H_*/
