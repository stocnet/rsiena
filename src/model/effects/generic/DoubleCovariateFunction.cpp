/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleCovariateFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * DoubleCovariateFunction.
 *****************************************************************************/

#include <stdexcept>
#include <math.h> /* round */
#include "DoubleCovariateFunction.h"
#include "data/Data.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/LongitudinalData.h"
#include "model/State.h"

using namespace std;

namespace siena
{

/**
 * Creates a new predicate.
 */
DoubleCovariateFunction::DoubleCovariateFunction(std::string covariateName1,
						std::string covariateName2) : AlterFunction()
{
	this->lFirstCovariateName = covariateName1;
	this->lSecondCovariateName = covariateName2;
	this->lpFirstConstantCovariate = 0;
	this->lpSecondConstantCovariate = 0;
	this->lpFirstChangingCovariate = 0;
	this->lpSecondChangingCovariate = 0;
	this->lpFirstBehaviorData = 0;
	this->lpSecondBehaviorData = 0;
	this->lFirstValues = 0;
	this->lSecondValues = 0;
}


/**
 * Initializes this predicate.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DoubleCovariateFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lperiod = period;

	this->lpFirstConstantCovariate =
				pData->pConstantCovariate(this->lFirstCovariateName);
	this->lpFirstChangingCovariate =
				pData->pChangingCovariate(this->lFirstCovariateName);
	this->lpFirstBehaviorData =
				pData->pBehaviorData(this->lFirstCovariateName);
	this->lFirstValues = pState->behaviorValues(this->lFirstCovariateName);
// these are the non-centered values

	if (!this->lpFirstConstantCovariate &&
		!this->lpFirstChangingCovariate &&
		!(this->lpFirstBehaviorData && this->lFirstValues))
	{
		throw logic_error("Covariate or dependent behavior variable '" +
			this->lFirstCovariateName +
			"' expected for first variable.");
	}

	this->lpSecondConstantCovariate =
				pData->pConstantCovariate(this->lSecondCovariateName);
	this->lpSecondChangingCovariate =
				pData->pChangingCovariate(this->lSecondCovariateName);
	this->lpSecondBehaviorData =
				pData->pBehaviorData(this->lSecondCovariateName);
	this->lSecondValues = pState->behaviorValues(this->lSecondCovariateName);
// the non-centered values

	if (!this->lpSecondConstantCovariate &&
		!this->lpSecondChangingCovariate &&
		!(this->lpSecondBehaviorData && this->lSecondValues))
	{
		throw logic_error("Covariate or dependent behavior variable '" +
			this->lSecondCovariateName +
			"' expected for second variable.");
	}
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling DoubleCovariateFunction::value(...).
 */
void DoubleCovariateFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
}


/**
 * Returns the first covariate value for the given actor.
 * For behavior, this is the non-centered value.
 */
double DoubleCovariateFunction::firstCovariateValue(int i) const
{
	double value = 0;

	if (this->lpFirstConstantCovariate)
	{
		value = this->lpFirstConstantCovariate->value(i);
	}
	else if (this->lpFirstChangingCovariate)
	{
		value = this->lpFirstChangingCovariate->value(i, this->lperiod);
	}
	else
	{
		value = this->lFirstValues[i];
	}
	return value;
}

/**
 * Returns the second covariate value for the given actor.
 * For behavior, this is the non-centered value.
 */
double DoubleCovariateFunction::secondCovariateValue(int i) const
{
	double value = 0;

	if (this->lpSecondConstantCovariate)
	{
		value = this->lpSecondConstantCovariate->value(i);
	}
	else if (this->lpSecondChangingCovariate)
	{
		value = this->lpSecondChangingCovariate->value(i, this->lperiod);
	}
	else
	{
		value = this->lSecondValues[i];
	}
	return value;
}


/**
 * Returns the first covariate number of cases.
 */
int DoubleCovariateFunction::firstCovariateN() const
{
	int ncov = 0;

	if (this->lpFirstConstantCovariate)
	{
		ncov = this->lpFirstConstantCovariate->covariateN();
	}
	else if (this->lpFirstChangingCovariate)
	{
		ncov = this->lpFirstChangingCovariate->covariateN();
	}
	else
	{
		ncov = this->lpFirstBehaviorData->observationCount();
	}
	return ncov;
}

/**
 * Returns the second covariate number of cases.
 */
int DoubleCovariateFunction::secondCovariateN() const
{
	int ncov = 0;

	if (this->lpSecondConstantCovariate)
	{
		ncov = this->lpSecondConstantCovariate->covariateN();
	}
	else if (this->lpSecondChangingCovariate)
	{
		ncov = this->lpSecondChangingCovariate->covariateN();
	}
	else
	{
		ncov = this->lpSecondBehaviorData->observationCount();
	}
	return ncov;
}

/**
 * Returns the first covariate value for the given actor, rounded to integer.
 * For behavior, this is the non-centered value.
 */
int DoubleCovariateFunction::firstCovariateIntValue(int i) const
{
	return int(round(this->firstCovariateValue(i)));
}

/**
 * Returns the second covariate value for the given actor, rounded to integer.
 * For behavior, this is the non-centered value.
 */
int DoubleCovariateFunction::secondCovariateIntValue(int i) const
{
	return int(round(this->secondCovariateValue(i)));
}


/**
 * Returns if the first covariate value for the given actor is missing.
 */
bool DoubleCovariateFunction::firstMissing(int i) const
{
	bool missing = false;

	if (this->lpFirstConstantCovariate)
	{
		missing = this->lpFirstConstantCovariate->missing(i);
	}
	else if (this->lpFirstChangingCovariate)
	{
		missing = this->lpFirstChangingCovariate->missing(i, this->lperiod);
	}
	else
	{
		missing = this->lpFirstBehaviorData->missing(this->lperiod, i);
	}

	return missing;
}

/**
 * Returns if the second covariate value for the given actor is missing.
 */
bool DoubleCovariateFunction::secondMissing(int i) const
{
	bool missing = false;

	if (this->lpSecondConstantCovariate)
	{
		missing = this->lpSecondConstantCovariate->missing(i);
	}
	else if (this->lpSecondChangingCovariate)
	{
		missing = this->lpSecondChangingCovariate->missing(i, this->lperiod);
	}
	else
	{
		missing = this->lpSecondBehaviorData->missing(this->lperiod, i);
	}

	return missing;
}

/**
 * Returns the first covariate minimum value.
 * For behavior, this is the minimum non-centered value.
 */
double DoubleCovariateFunction::firstCovariateMinimum() const
{
	double mini = 0;

	if (this->lpFirstConstantCovariate)
	{
		mini = this->lpFirstConstantCovariate->min();
	}
	else if (this->lpFirstChangingCovariate)
	{
		mini = this->lpFirstChangingCovariate->min();
	}
	else
	{
		mini = this->lpFirstBehaviorData->min();
	}
	return mini;
}


/**
 * Returns the first covariate maximum value.
 * For behavior, this is the maximum non-centered value.
 */
double DoubleCovariateFunction::firstCovariateMaximum() const
{
	double maxi = 0;

	if (this->lpFirstConstantCovariate)
	{
		maxi = this->lpFirstConstantCovariate->max();
	}
	else if (this->lpFirstChangingCovariate)
	{
		maxi = this->lpFirstChangingCovariate->max();
	}
	else
	{
		maxi = this->lpFirstBehaviorData->max();
	}
	return maxi;
}


/**
 * Returns the second covariate minimum value.
 * For behavior, this is the minimum non-centered value.
 */
double DoubleCovariateFunction::secondCovariateMinimum() const
{
	double mini = 0;

	if (this->lpSecondConstantCovariate)
	{
		mini = this->lpSecondConstantCovariate->min();
	}
	else if (this->lpSecondChangingCovariate)
	{
		mini = this->lpSecondChangingCovariate->min();
	}
	else
	{
		mini = this->lpSecondBehaviorData->min();
	}
	return mini;
}

/**
 * Returns the second covariate maximum value.
 * For behavior, this is the maximum non-centered value.
 */
double DoubleCovariateFunction::secondCovariateMaximum() const
{
	double maxi = 0;

	if (this->lpSecondConstantCovariate)
	{
		maxi = this->lpSecondConstantCovariate->max();
	}
	else if (this->lpSecondChangingCovariate)
	{
		maxi = this->lpSecondChangingCovariate->max();
	}
	else
	{
		maxi = this->lpSecondBehaviorData->max();
	}
	return maxi;
}

}
