/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreeActivitySqrtEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeActivitySqrtEffect class.
 *****************************************************************************/

#ifndef OUTDEGREEACTIVITYSQRTEFFECT_H_
#define OUTDEGREEACTIVITYSQRTEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;
class Network;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * This class defines the outdegree activity (sqrt) effect defined by
 * s_i(x)= x_{i+}^1.5. The corresponding statistic is sum_i x_{i+}^1.5.
 */
class OutdegreeActivitySqrtEffect : public NetworkEffect
{
friend class BothDegreesEffect;

public:
	OutdegreeActivitySqrtEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	virtual double endowmentStatistic(Network * pLostTieNetwork);

private:
	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*OUTDEGREEACTIVITYSQRTEFFECT_H_*/
