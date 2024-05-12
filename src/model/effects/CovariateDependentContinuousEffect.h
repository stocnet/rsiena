/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDependentContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDependentContinuousEffect class.
 *****************************************************************************/

#ifndef COVARIATEDEPENDENTCONTINUOUSEFFECT_H_
#define COVARIATEDEPENDENTCONTINUOUSEFFECT_H_

#include "ContinuousEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;
class ChangingCovariate;
class ContinuousVariable;
class BehaviorLongitudinalData;
class ContinuousLongitudinalData;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for continuous behavior effects depending on an individual
 * covariate (constant, changing, or other dependent behavior variable).
 */
class CovariateDependentContinuousEffect : public ContinuousEffect
{
public:
	CovariateDependentContinuousEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	double covariateValue(int i) const;
	bool missingCovariate(int i, int observation) const;
	bool missingCovariateEitherEnd(int i, int observation) const;

private:
	ConstantCovariate * lpConstantCovariate;
	ChangingCovariate * lpChangingCovariate;
	BehaviorLongitudinalData * lpBehaviorData;
	ContinuousLongitudinalData * lpContinuousData;
	
	// The current value of a (discrete or continuous) interacting behavior 
	// variable per each actor.
	// This array is 0 for covariate-based effects.
	
	const int * lvalues {};
	const double * lcontinuousValues {};
};

}

#endif /*COVARIATEDEPENDENTCONTINUOUSEFFECT_H_*/
