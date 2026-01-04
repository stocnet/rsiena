/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: alterscovariateavrecalteffect.h
 *
 * Description: This file contains the definition of the
 * AltersCovariateAvRecAltEffect class.
 * Made on the basis of AltersCovariateAvAltEffect.
 *****************************************************************************/

#ifndef ALTERSCOVARIATEAVRECALTEFFECT_H_
#define ALTERSCOVARIATEAVRECALTEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Alters covariate average alter effect 
 */
class AltersCovariateAvRecAltEffect :
public CovariateAndNetworkBehaviorEffect
{
public:
	AltersCovariateAvRecAltEffect(const EffectInfo * pEffectInfo, 
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

#endif /*ALTERSCOVARIATEAVRECALTEFFECT_H_*/
