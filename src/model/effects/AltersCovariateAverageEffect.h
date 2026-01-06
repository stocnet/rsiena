/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAverageEffect.h
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
 * Alters covariate average alter effect 
 */
class AltersCovariateAverageEffect :
public CovariateAndNetworkBehaviorEffect
{
public:
	AltersCovariateAverageEffect(const EffectInfo * pEffectInfo, 
			bool divide, bool same, bool outgoing);

	virtual void preprocessEgo(int ego);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

protected:

private:
	// divide indicates whether there will be division by the outdegree
	bool ldivide{};
	// same indicates whether the influence will be from same-X alters
	bool lsame{};
	// outgoing indicates whether the influence will be from out-alters
	bool loutgoing{};
	double lTotalAlterValue{};
};

}

#endif /*AltersCovariateAverageEffect_H_*/
