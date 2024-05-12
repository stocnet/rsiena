/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WXXClosureEffect.h
 *
 * Description: This file contains the definition of the
 * WXXClosureEffect class.
 *****************************************************************************/

#ifndef WXXCLOSUREEFFECT_H_
#define WXXCLOSUREEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * WX => X closure of covariate effect (see manual).
 */
class WXXClosureEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	WXXClosureEffect(const EffectInfo * pEffectInfo);
	virtual ~WXXClosureEffect();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	void calculateSums(int i, const Network * pNetwork, double * sums) const;

	// For a fixed i, this variable stores the value of sum_h w_{ih} x_{hj} for
	// each j.
	double * lsums {};
};

}

#endif /*WXXCLOSUREEFFECT_H_*/
