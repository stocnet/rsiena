/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WWXClosureEffect.h
 *
 * Description: This file contains the definition of the
 * WWXClosureEffect class.
 *****************************************************************************/

#ifndef WWXCLOSUREEFFECT_H_
#define WWXCLOSUREEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * WW => X closure of covariate effect (see manual);
 * different values of directionality of W tie out1, out2
 */
class WWXClosureEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	WWXClosureEffect(const EffectInfo * pEffectInfo, bool out1, bool out2);
	virtual ~WWXClosureEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	void calculateSums(int i, const Network * pNetwork, double * sums) const;

	// For a fixed i, this variable stores the value for each j of
	// out1 = out2 = TRUE:  sum_h w_{ih} w_{hj}
	// out1 = out2 = FALSE: sum_h w_{hi} w_{jh}
	// out1 = FALSE, out2 = TRUE: sum_h w_{hi} w_{hj}

	double * lsums {};
	bool lout1 {};
	bool lout2 {};

};

}

#endif /*WWXCLOSUREEFFECT_H_*/
