/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InAltersCovariateAverageBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * InAltersCovariateAverageEffect class.
 *****************************************************************************/

#ifndef INALTERSCOVARIATEAVERAGEEFFECT_H_
#define INALTERSCOVARIATEAVERAGEEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * In-alters covariate average effect (see manual)
 *
 */
class InAltersCovariateAverageEffect :
public CovariateAndNetworkBehaviorEffect
{
public:
	InAltersCovariateAverageEffect(const EffectInfo * pEffectInfo, bool divide);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

protected:

private:
	// divide indicates whether there will be division by the outdegree
	bool ldivide {};
};

}

#endif /*INALTERSCOVARIATEAVERAGEEFFECT_H_*/
