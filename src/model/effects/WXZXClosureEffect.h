/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WXZXClosureEffect.h
 *
 * Description: Definition of the WXZ => X multilevel four-cycle closure
 * effect for a bipartite dependent network X, a first-mode dyadic covariate
 * W, and a second-mode dyadic covariate Z.
 *****************************************************************************/

#ifndef WXZXCLOSUREEFFECT_H_
#define WXZXCLOSUREEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

class ConstantDyadicCovariate;
class ChangingDyadicCovariate;
class DyadicCovariateValueIterator;

/**
 * WXZ => X multilevel four-cycle closure.
 *
 * For a focal first-mode actor i and second-mode alter j, the contribution is
 * sum_{h,k} w_ih x_hk z_kj.
 */
class WXZXClosureEffect : public NetworkEffect
{
public:
	WXZXClosureEffect(const EffectInfo * pEffectInfo);
	virtual ~WXZXClosureEffect();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	virtual void initializeStatisticCalculation();
	virtual void cleanupStatisticCalculation();

private:
	DyadicCovariateValueIterator wRowValues(int i) const;
	DyadicCovariateValueIterator zRowValues(int k) const;

	ConstantDyadicCovariate * lpWConstant {};
	ChangingDyadicCovariate * lpWChanging {};
	ConstantDyadicCovariate * lpZConstant {};
	ChangingDyadicCovariate * lpZChanging {};

	// For fixed ego i: sum_{h,k} w_ih x_hk z_kj, for every second-mode j.
	double * lsums {};
	bool lexcludeMissings {};
};

}

#endif /* WXZXCLOSUREEFFECT_H_ */
