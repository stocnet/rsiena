/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameInCovariateFourCyclesEffect.h
 *
 * Description: This file contains the definition of the
 * SameInCovariateFourCyclesEffect class.
 *****************************************************************************/

#ifndef SAMEINCOVARIATEFOURCYCLESEFFECT_H_
#define SAMEINCOVARIATEFOURCYCLESEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;


// ----------------------------------------------------------------------------
// Section: SameInCovariateFourCyclesEffect class
// ----------------------------------------------------------------------------

/**
 * Same and different covariate FourCycles effects (see manual).
 */
class SameInCovariateFourCyclesEffect : public CovariateDependentNetworkEffect
{
public:
	SameInCovariateFourCyclesEffect(const EffectInfo * pEffectInfo);
	virtual ~SameInCovariateFourCyclesEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego, const Network * pSummationTieNetwork);

private:
	// Indicates if the effect is used for a two-mode or one-mode network
	bool lTwoMode {};

	void countThreePaths(int i, const Network * pNetwork, long int * counters, 
			bool dropMissings) const;
	// For a fixed i, this variable stores the number of three-paths
	// i -> h <- k -> j per each j.

	long int * lcounters {};

	// Indicates if the square root of the number of four-cycles has to
	// be taken.
	bool lroot {};
	// Lookup table for fast square root calculations:
	SqrtTable * lpSqrtTable;

	// The number of 4-cycles the ego is currently involved in.
	int lcurrentCycleCount {};
};

}

#endif /*SAMEINCOVARIATEFOURCYCLESEFFECT_H_*/
