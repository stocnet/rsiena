/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateGwdspEffect.h
 *
 * Description: This file contains the declaration of the class
 * CovariateGwdspEffect.
 *
 * This is the covariate/behaviour-weighted geometrically weighted dyadwise
 * shared partners network (selection) effect -- the network-equation
 * building block whose behaviour-equation twin is TotalGwdspAlterEffect.
 * Per potential tie i -> j its change contribution is
 *   sum_{h in in/outTies(j), h != i} z_h * ( w(S_ih + 1) - w(S_ih) ),
 * where S_ih = # shared partners of i and h (dyadwise), w is the GW
 * transform, and z_h is the covariate/behaviour value of h. The focal
 * alter j is not weighted directly: adding tie i -> j raises S_ih for
 * exactly the h that also connect to j, so it is those h's behaviour that
 * enters (see docs/design). Multiply by egoX_nc (z_i) via a user
 * interaction to obtain the full z_i * sum_h z_h * w(S_ih) statistic.
 *****************************************************************************/

#ifndef COVARIATEGWDSPEFFECT_H_
#define COVARIATEGWDSPEFFECT_H_

#include "model/effects/CovariateDependentNetworkEffect.h"
#include <vector>

namespace siena
{

class EgocentricConfigurationTable;
class NetworkCache;

/**
 * Covariate-weighted geometrically weighted dyadwise shared partners effect.
 */
class CovariateGwdspEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateGwdspEffect(const EffectInfo * pEffectInfo,
			EgocentricConfigurationTable * (NetworkCache::*pTable)() const,
			bool forward, bool nc);
	virtual void initialize(const Data * pData, State * pState, int period,
			Cache * pCache);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);

private:

	NetworkCache * lpNetworkCache;
	EgocentricConfigurationTable * (NetworkCache::*lpTable)() const;
	double lparameter {};
	std::vector<double> lcumulativeWeight;
	bool lforward {};
	bool lnc {};
	double lweight {};
	double lexpmweight {};
	double lexpfactor {};
	EgocentricConfigurationTable * lpInitialisedTable;
};

}

#endif /*COVARIATEGWDSPEFFECT_H_*/
