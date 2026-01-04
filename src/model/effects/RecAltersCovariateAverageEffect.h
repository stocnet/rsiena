/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: RecAltersCovariateAverageEffect.h
 *
 * Description: This file contains the definition of the
 * RecAltersCovariateAverageEffect class.
 * Made on the basis of AltersCovariateAvAltEffect.
 *****************************************************************************/

#ifndef RECALTERSCOVARIATEAVERAGEEFFECT_H_
#define RECALTERSCOVARIATEAVERAGEEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Alters covariate average alter effect 
 */
class RecAltersCovariateAverageEffect :
public CovariateAndNetworkBehaviorEffect
{
public:
	RecAltersCovariateAverageEffect(const EffectInfo * pEffectInfo, 
			bool divide, bool same);
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
	double lTotalAlterValue{};
};

}

#endif /*RECALTERSCOVARIATEAVERAGEEFFECT_H_*/
