/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: XXWClosureEffect.h
 *
 * Description: This file contains the definition of the
 * XXWClosureEffect class.
 *****************************************************************************/

#ifndef XXWCLOSUREEFFECT_H_
#define XXWCLOSUREEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * XW => X closure of covariate effect (see manual).
 */
class XXWClosureEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	XXWClosureEffect(const EffectInfo * pEffectInfo, bool outst, bool inst);
	virtual ~XXWClosureEffect();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	void calculateOutStarSums(int i, const Network * pNetwork,
		double * sums) const;
	void calculateInStarSums(int i, const Network * pNetwork,
		double * sums) const;
// note that this function differs from XWXClosureEffect.calculateInStarSums

	// The network this effect is associated with
	const Network * lpNetwork;

	// For a fixed i, this variable stores the value of sum_h x_{hi} w_{hj} for
	// each j.
	double * loutStarSums {};

	// For a fixed i, this variable stores the value of sum_h w_{ih} x_{jh} for
	// each j.
	double * linStarSums {};

	// divide indicates whether there will be division by the indegree
	bool loutst {}; // contribution from outstars
	bool linst {}; // contribution from instars
};

}

#endif /*XXWCLOSUREEFFECT_H_*/
