/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: FourCyclesEffect.h
 *
 * Description: This file contains the definition of the
 * FourCyclesEffect class.
 *****************************************************************************/

#ifndef FOURCYCLESEFFECT_H_
#define FOURCYCLESEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;


// ----------------------------------------------------------------------------
// Section: FourCyclesEffect class
// ----------------------------------------------------------------------------

/**
 * 4-cycles effect (see manual).
 */
class FourCyclesEffect : public NetworkEffect
{
public:
	FourCyclesEffect(const EffectInfo * pEffectInfo, bool TwoMode);
	virtual ~FourCyclesEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the effect is used for a two-mode or one-mode network
	bool lTwoMode {};

	void countThreePaths(int i, const Network * pNetwork, long int * counters)
		const;

	// For a fixed i, this variable stores the number of three-paths
	// i -> h <- k -> j per each j.

	long int * lcounters {};

	// Indicates if the square root of the number of four-cycles has to
	// be taken.

	bool lroot {};

	// Lookup table for fast square root calculations
	SqrtTable * lpSqrtTable;

	// The number of 4-cycles the ego is currently involved in.
	int lcurrentCycleCount {};
};

}

#endif /*FOURCYCLESEFFECT_H_*/
