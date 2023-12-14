/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DenseTriadsEffect.h
 *
 * Description: This file contains the declaration of the class
 * DenseTriadsEffect.
 *****************************************************************************/

#ifndef DENSETRIADSEFFECT_H_
#define DENSETRIADSEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the dense triads effect. See the manual for effect
 * definitions.
 */
class DenseTriadsEffect : public NetworkEffect
{
public:
	DenseTriadsEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual void initializeStatisticCalculation();
	virtual void onNextEgo(int ego);
	virtual double tieStatistic(int alter);
	virtual void cleanupStatisticCalculation();

private:
	int ldensity {};

	// A helper array of marks for statistic calculations
	int * lmark {};

	// Given an ego i
	// mark[h] = baseMark + 2 if there are mutual ties between i and h,
	// mark[h] = baseMark + 1 if only one of the mutual ties is present,
	// mark[h] <= baseMark otherwise.

	int lbaseMark {};
};

}

#endif /*DENSETRIADSEFFECT_H_*/
