/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * EgoFunction.
 *****************************************************************************/

#include "EgoFunction.h"

namespace siena
{

/**
 * Creates a new function which uses ego's value
 */
EgoFunction::EgoFunction(AlterFunction * pFirstFunction)
{
	this->lpFirstFunction = pFirstFunction;
}


/**
 * Deallocates this function.
 */
EgoFunction::~EgoFunction()
{
	delete this->lpFirstFunction;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void EgoFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpFirstFunction->initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling DifferenceFunction::value(...).
 */
void EgoFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	this->lpFirstFunction->preprocessEgo(ego);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 * Note that this->ego() is used instead of alter.
 */
double EgoFunction::value(int alter) const
{
	return this->lpFirstFunction->value(this->ego());
}

}
