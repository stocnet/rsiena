/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDependentBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDependentBehaviorEffect class.
 *****************************************************************************/

#ifndef COVARIATEDEPENDENTBEHAVIOREFFECT_H_
#define COVARIATEDEPENDENTBEHAVIOREFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;
class ChangingCovariate;
class BehaviorVariable;
class BehaviorLongitudinalData;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for behavior effects depending on an individual
 * covariate (constant, changing, or other dependent behavior variable).
 */
class CovariateDependentBehaviorEffect : public BehaviorEffect
{
public:
	CovariateDependentBehaviorEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	double covariateMean() const;
	double covariateValue(int i) const;
	bool missingCovariate(int i, int observation) const;
	bool missingCovariateEitherEnd(int i, int observation) const;

private:
	ConstantCovariate * lpConstantCovariate;
	ChangingCovariate * lpChangingCovariate;
	BehaviorLongitudinalData * lpBehaviorData;

	// The current value of an interacting behavior variable
	// per each actor.
	// This array is 0 for covariate-based effects.

	const int * linteractionValues {};
};

}

#endif /*COVARIATEDEPENDENTBEHAVIOREFFECT_H_*/
