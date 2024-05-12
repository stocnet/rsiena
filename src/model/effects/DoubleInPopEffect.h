/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleInPopEffect.h
 *
 * Description: This file contains the definition of the
 * DoubleInPopEffect class.
 *****************************************************************************/

#ifndef DOUBLEINPOPEFFECT_H_
#define DOUBLEINPOPEFFECT_H_

#include "MixedNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: DoubleInPopEffect class
// ----------------------------------------------------------------------------

/**
 * This class defines the double indegree popularity effect defined by
 * s_i(x)= w{ij} * (sum_h x_{hj} w{{hj}).
 */
class DoubleInPopEffect : public MixedNetworkEffect
{
public:
	DoubleInPopEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool lroot {};
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /*DOUBLEINPOPEFFECT_H_*/
