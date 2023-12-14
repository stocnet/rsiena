/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutOutDegreeAssortativityEffect.h
 *
 * Description: This file contains the definition of the
 * OutOutDegreeAssortativityEffect class.
 *****************************************************************************/

#ifndef OUTOUTDEGREEASSORTATIVITYEFFECT_H_
#define OUTOUTDEGREEASSORTATIVITYEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The out-out degree assortativity effect (see manual).
 */
class OutOutDegreeAssortativityEffect : public NetworkEffect
{
public:
	OutOutDegreeAssortativityEffect(const EffectInfo * pEffectInfo);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the square root of degrees must be used
	bool lroot {};

	// The sum of out-degrees (or their square roots) over all
	// out-neighbors of the ego

	double lneighborDegreeSum {};

	// The current out-degree of the ego
	int ldegree {};

	// The square root of the current out-degree of the ego
	double lsqrtDegree {};

	// sqrt(current out-degree of ego + 1)
	double lsqrtDegreePlus {};

	// sqrt(current out-degree of ego - 1)
	double lsqrtDegreeMinus {};

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*OUTOUTDEGREEASSORTATIVITYEFFECT_H_*/
