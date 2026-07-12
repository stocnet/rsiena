/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ActivityAlterEffect.h
 *
 * Description: This file contains the definition of the
 * ActivityAlterEffect class.
 *****************************************************************************/

#ifndef ACTIVITYALTEREFFECT_H_
#define ACTIVITYALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average and total Activity alter effect: ego's behavior weighted by the
 * (average) out-degree of its alters. This is the unweighted FF structural
 * control in behavior dynamics, analogous to popAlt (FB control).
 */
class ActivityAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	ActivityAlterEffect(const EffectInfo * pEffectInfo, bool divide, bool nc);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	double averageOutDegree(int i) const;
	// divide indicates whether there will be division by the outdegree of ego
	bool ldivide {};
	bool lnc {};
};

}

#endif /* ACTIVITYALTEREFFECT_H_ */
