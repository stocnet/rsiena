/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariatePredicate.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariatePredicate.
 *****************************************************************************/

#include <stdexcept>
#include "CovariatePredicate.h"
#include "data/Data.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"

using namespace std;

namespace siena
{

/**
 * Creates a new predicate.
 */
CovariatePredicate::CovariatePredicate(string covariateName) :
	NamedObject(covariateName)
{
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
	this->lpBehaviorData = 0;
	this->lvalues = 0;
}


/**
 * Initializes this predicate.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariatePredicate::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterPredicate::initialize(pData, pState, period, pCache);
	string name = this->name();

	this->lpConstantCovariate = pData->pConstantCovariate(name);
	this->lpChangingCovariate = pData->pChangingCovariate(name);
	this->lpBehaviorData = pData->pBehaviorData(name);
	this->lvalues = pState->behaviorValues(name);

	if (!this->lpConstantCovariate &&
		!this->lpChangingCovariate &&
		!(this->lpBehaviorData && this->lvalues))
	{
		throw logic_error("Covariate or dependent behavior variable '" +
			name +
			"' expected.");
	}
}


/**
 * Returns the covariate value for the given actor.
 */
double CovariatePredicate::covariateValue(int i) const
{
	double value = 0;

	if (this->lpConstantCovariate)
	{
		value = this->lpConstantCovariate->value(i);
	}
	else if (this->lpChangingCovariate)
	{
		value = this->lpChangingCovariate->value(i, this->period());
	}
	else
	{
		value = this->lvalues[i] - this->lpBehaviorData->overallMean();
	}

	return value;
}


/**
 * Returns if the covariate value for the given actor is missing.
 */
bool CovariatePredicate::missing(int i) const
{
	bool missing = false;

	if (this->lpConstantCovariate)
	{
		missing = this->lpConstantCovariate->missing(i);
	}
	else if (this->lpChangingCovariate)
	{
		missing = this->lpChangingCovariate->missing(i, this->period());
	}
	else
	{
		missing = this->lpBehaviorData->missing(this->period(), i);
	}

	return missing;
}

}
