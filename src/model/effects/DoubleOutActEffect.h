/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleOutActEffect.h
 *
 * Description: This file contains the definition of the
 * DoubleOutActEffect class.
 *****************************************************************************/

#ifndef DOUBLEOUTACTEFFECT_H_
#define DOUBLEOUTACTEFFECT_H_

#include "MixedNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: DoubleOutActEffect class
// ----------------------------------------------------------------------------

/**
 * This class defines the double outdegree activity effect defined by
 * s_i(x)= (sum_h x_{ih} w{{ih})^2.
 */
class DoubleOutActEffect : public MixedNetworkEffect
{
public:
	DoubleOutActEffect(const EffectInfo * pEffectInfo);

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

#endif /*DOUBLEOUTACTEFFECT_H_*/
