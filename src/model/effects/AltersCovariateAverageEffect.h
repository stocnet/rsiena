/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAverageBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * AltersCovariateAverageEffect class.
 *****************************************************************************/

#ifndef ALTERSCOVARIATEAVERAGEEFFECT_H_
#define ALTERSCOVARIATEAVERAGEEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Alters covariate average effect (see manual)
 *
 */
class AltersCovariateAverageEffect :
public CovariateAndNetworkBehaviorEffect
{
public:
	AltersCovariateAverageEffect(const EffectInfo * pEffectInfo, bool divide);

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

#endif /*ALTERSCOVARIATEAVERAGEEFFECT_H_*/
