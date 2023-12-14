/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreeContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeEffect class.
 *****************************************************************************/

#ifndef INDEGREECONTINUOUSEFFECT_H_
#define INDEGREECONTINUOUSEFFECT_H_

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
 * Indegree effect defined as the number of its inward neighbors (with 
 * respect to a certain network).
 */
class IndegreeContinuousEffect : public NetworkDependentContinuousEffect
{
public:
	IndegreeContinuousEffect(const EffectInfo * pEffectInfo, bool root);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);
	
private:
	// Indicates if the square root of indegrees must be used
	bool lroot {};

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;	
};

}

#endif /*INDEGREECONTINUOUSEFFECT_H_*/

