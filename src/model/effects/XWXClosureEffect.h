/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: XWXClosureEffect.h
 *
 * Description: This file contains the definition of the
 * XWXClosureEffect class.
 *****************************************************************************/

#ifndef XWXCLOSUREEFFECT_H_
#define XWXCLOSUREEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * XW => X closure of covariate effect (see manual).
 */
class XWXClosureEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	XWXClosureEffect(const EffectInfo * pEffectInfo, bool tp, bool inst);
	virtual ~XWXClosureEffect();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	void calculateTwoPathSums(int i, const Network * pNetwork,
		double * sums) const;
	void calculateInStarSums(int i, const Network * pNetwork,
		double * sums) const;
// note that this function differs from XXWClosureEffect.calculateInStarSums

	// For a fixed i, this variable stores the value of sum_h x_{ih} w_{hj} for
	// each j.
	double * ltwoPathSums {};

	// For a fixed i, this variable stores the value of sum_h x_{ih} w_{jh} for
	// each j.
	double * linStarSums {};

	// divide indicates whether there will be division by the indegree
	bool ltp {}; // contribution from twopaths
	bool linst {}; // contribution from instars
};

}

#endif /*XWXCLOSUREEFFECT_H_*/
