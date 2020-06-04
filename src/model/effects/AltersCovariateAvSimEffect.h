/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAvSimEffect.h
 *
 * Description: This file contains the definition of the
 * AltersCovariateAvSimEffect class.
 *****************************************************************************/

#ifndef ALTERSCOVARIATEAVSIMEFFECT_H_
#define ALTERSCOVARIATEAVSIMEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Alters covariate average similarity effect (not in manual)
 */
class AltersCovariateAvSimEffect : 
public CovariateAndNetworkBehaviorEffect
{
public:
	AltersCovariateAvSimEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

protected:

private:
};

}

#endif /*ALTERSCOVARIATEAVSIMEFFECT_H_*/
