/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConditionalFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * ConditionalFunction.
 *****************************************************************************/

#include "ConditionalFunction.h"
#include "model/effects/generic/AlterPredicate.h"

namespace siena
{

/**
 * Creates a new conditional function.
 * @param[in] pPredicate the predicate to be tested
 * @param[in] pIfFunction the function to be called if the predicate holds
 * @param[in] pElseFunction the function to be called if the predicate does not
 * hold
 */
ConditionalFunction::ConditionalFunction(AlterPredicate * pPredicate,
	AlterFunction * pIfFunction,
	AlterFunction * pElseFunction)
{
	this->lpPredicate = pPredicate;
	this->lpIfFunction = pIfFunction;
	this->lpElseFunction = pElseFunction;
}


/**
 * Deallocates this function.
 */
ConditionalFunction::~ConditionalFunction()
{
	delete this->lpPredicate;
	delete this->lpIfFunction;
	delete this->lpElseFunction;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ConditionalFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpPredicate->initialize(pData, pState, period, pCache);

	if (this->lpIfFunction)
	{
		this->lpIfFunction->initialize(pData, pState, period, pCache);
	}

	if (this->lpElseFunction)
	{
		this->lpElseFunction->initialize(pData, pState, period, pCache);
	}
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling ConditionalFunction::value(...).
 */
void ConditionalFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	this->lpPredicate->preprocessEgo(ego);

	if (this->lpIfFunction)
	{
		this->lpIfFunction->preprocessEgo(ego);
	}

	if (this->lpElseFunction)
	{
		this->lpElseFunction->preprocessEgo(ego);
	}
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double ConditionalFunction::value(int alter)
{
	double value = 0;

	if (this->lpPredicate->value(alter))
	{
		if (this->lpIfFunction)
		{
			value = this->lpIfFunction->value(alter);
		}
	}
	else if (this->lpElseFunction)
	{
		value = this->lpElseFunction->value(alter);
	}

	return value;
}

}
