/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreeContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeEffect class.
 *****************************************************************************/

#ifndef OUTDEGREECONTINUOUSEFFECT_H_
#define OUTDEGREECONTINUOUSEFFECT_H_

#include "NetworkDependentContinuousEffect.h"

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
 * Outdegree effect defined as the number of its outward neighbors (with 
 * respect to a certain network).
 */
class OutdegreeContinuousEffect : public NetworkDependentContinuousEffect
{
public:
	OutdegreeContinuousEffect(const EffectInfo * pEffectInfo, bool root);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);
	
private:
	// Indicates if the square root of indegrees must be used
	bool lroot {};

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;	
};

}

#endif /*OUTDEGREECONTINUOUSEFFECT_H_*/

