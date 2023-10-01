/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocalFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * ReciprocalFunction.
 *****************************************************************************/

#include "ReciprocalFunction.h"

namespace siena
{

/**
 * Creates a new function as the reciprocal of the given function.
 */
ReciprocalFunction::ReciprocalFunction(AlterFunction * pFunction)
{
	this->lpFunction = pFunction;
}


/**
 * Deallocates this function.
 */
ReciprocalFunction::~ReciprocalFunction()
{
	delete this->lpFunction;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ReciprocalFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpFunction->initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego.
 */
void ReciprocalFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	this->lpFunction->preprocessEgo(ego);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double ReciprocalFunction::value(int alter) const
{
	double fvalue = this->lpFunction->value(alter);
	if (!(fvalue == 0))
	{
		fvalue = 1/fvalue;
	}
	return fvalue;
}

}
