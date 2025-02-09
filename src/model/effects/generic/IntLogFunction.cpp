/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: https://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IntLogFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * IntLogFunction.
 *****************************************************************************/

#include "IntLogFunction.h"
#include "utils/LogTable.h"

namespace siena
{

IntLogFunction::IntLogFunction(AlterFunction * pFunction)
{
	this->lpFunction = pFunction;
	this->lpLogTable = LogTable::instance();
}


/**
 * Deallocates this function.
 */
IntLogFunction::~IntLogFunction()
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
void IntLogFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpFunction->initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling IntLogFunction::value(...).
 */
void IntLogFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	this->lpFunction->preprocessEgo(ego);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double IntLogFunction::value(int alter) const
{
	return this->lpLogTable->log(this->lpFunction->value(alter));
}

}
