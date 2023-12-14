/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DistanceTwoEffect.h
 *
 * Description: This file contains the declaration of the class
 * DistanceTwoEffect.
 *****************************************************************************/

#ifndef DISTANCETWOEFFECT_H_
#define DISTANCETWOEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines a generalization of the "number of distances two" effect,
 * where one can specify the number of two-paths required between a pair of
 * actors to qualify for a distance two pair. The manual defines the effect for
 * parameter values 1 and 2.
 */
class DistanceTwoEffect : public NetworkEffect
{
public:
	DistanceTwoEffect(const EffectInfo * pEffectInfo,
		int requiredTwoPathCount);

	virtual double calculateContribution(int alter) const;
	virtual double endowmentStatistic(Network * pLostTieNetwork);

protected:
	virtual void initializeStatisticCalculation();
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);
	virtual void cleanupStatisticCalculation();

private:
	// A helper array of marks for statistic calculation
	int * lmark {};
	int lrequiredTwoPathCount {};
};

}

#endif /*DISTANCETWOEFFECT_H_*/
