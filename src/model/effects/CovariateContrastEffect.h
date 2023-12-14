/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateContrastEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateContrastEffect class.
 *****************************************************************************/

#ifndef COVARIATECONTRASTEFFECT_H_
#define COVARIATECONTRASTEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{

/**
 * Covariate Contrast effect defined as the product of ego's behavior with the difference
 * (absolute / positive part / negative part) between covariate and outdegree.
 */
 
class CovariateContrastEffect : public CovariateAndNetworkBehaviorEffect
{
public:
	CovariateContrastEffect(const EffectInfo * pEffectInfo, 
			const bool plus, const bool minus);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);

private:
	bool lplus {};
	bool lminus {};
};

}

#endif /* COVARIATECONTRASTEFFECT_H_*/
