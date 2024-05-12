/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageGroupEgoEffect.h
 *
 * Description: This file contains the definition of the
 * AverageGroupEgoEffect class.
 *****************************************************************************/

#ifndef AVERAGEGROUPEGOEFFECT_H_
#define AVERAGEGROUPEGOEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: AverageGroupEgoEffect class
// ----------------------------------------------------------------------------


/**
 * Group-average ego effect.
 */
class AverageGroupEgoEffect : public CovariateDependentNetworkEffect
{
public:
	AverageGroupEgoEffect(const EffectInfo * pEffectInfo);

	void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;
	virtual bool egoEffect() const;


protected:
	virtual double tieStatistic(int alter);
	virtual double egoStatistic(int ego, const Network * pNetwork);

private:
	// lcentermean = whether to center about the general mean
//	bool lcenterMean {};
	// if not lcenter, centering is about the following value
//	double lcenteringValue {};
	double lnm {};
	double lGroupMean {};
//	double loverallCenterMean {};
	BehaviorLongitudinalData * lpBehaviorData {};
	int lperiod {};
};

}

#endif /*AVERAGEGROUPEGOEFFECT_H_*/
