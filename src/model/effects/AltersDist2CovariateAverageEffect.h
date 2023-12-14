/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersDist2CovariateAverageBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * AltersDist2CovariateAverageEffect class.
 *****************************************************************************/

#ifndef ALTERSDIST2COVARIATEAVERAGEEFFECT_H_
#define ALTERSDIST2COVARIATEAVERAGEEFFECT_H_

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
class AltersDist2CovariateAverageEffect :
public CovariateAndNetworkBehaviorEffect
{
public:
	AltersDist2CovariateAverageEffect(const EffectInfo * pEffectInfo, bool divide1, bool divide2);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

protected:

private:
	bool ldivide1{};
	// Indicates whether there will be division by the outdegree of ego
	bool ldivide2{};
	// Indicates whether there will be division by the outdegree of alter
};

}

#endif /*ALTERSDIST2COVARIATEAVERAGEEFFECT_H_*/
