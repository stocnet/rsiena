/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateTotSimEffect.h
 *
 * Description: This file contains the definition of the
 * AltersCovariateTotSimEffect class.
 *****************************************************************************/

#ifndef ALTERSCOVARIATETOTSIMEFFECT_H_
#define ALTERSCOVARIATETOTSIMEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Alters covariate average similarity effect (not in manual)
 */
class AltersCovariateTotSimEffect : 
public CovariateAndNetworkBehaviorEffect
{
public:
	AltersCovariateTotSimEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

protected:

private:
};

}

#endif /*ALTERSCOVARIATETOTSIMEFFECT_H_*/
