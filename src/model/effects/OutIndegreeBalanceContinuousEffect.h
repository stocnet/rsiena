/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutIndegreeBalaceContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * OutIndegreeBalanceContinuousEffect class.
 *****************************************************************************/

#ifndef OUTINDEGREEBALANCECONTINUOUSEFFECT_H_
#define OUTINDEGREEBALANCECONTINUOUSEFFECT_H_

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
class OutIndegreeBalanceContinuousEffect : public NetworkDependentContinuousEffect
{
public:
	OutIndegreeBalanceContinuousEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor);
	virtual double egoStatistic(int ego, double * currentValues);
	
private:
    double outIndegreeBalance(int ego) const;
};

}

#endif /*OUTINDEGREEBALANCECONTINUOUSEFFECT_H_*/

