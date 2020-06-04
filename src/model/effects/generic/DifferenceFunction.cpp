/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DifferenceFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * DifferenceFunction.
 *****************************************************************************/

#include "DifferenceFunction.h"

namespace siena
{

/**
 * Creates a new function as a difference between the value of the two
 * give functions.
 */
DifferenceFunction::DifferenceFunction(AlterFunction * pFirstFunction,
	AlterFunction * pSecondFunction)
{
	this->lpFirstFunction = pFirstFunction;
	this->lpSecondFunction = pSecondFunction;
}


/**
 * Deallocates this function.
 */
DifferenceFunction::~DifferenceFunction()
{
	delete this->lpFirstFunction;
	delete this->lpSecondFunction;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DifferenceFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpFirstFunction->initialize(pData, pState, period, pCache);
	this->lpSecondFunction->initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling DifferenceFunction::value(...).
 */
void DifferenceFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	this->lpFirstFunction->preprocessEgo(ego);
	this->lpSecondFunction->preprocessEgo(ego);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double DifferenceFunction::value(int alter)
{
	return this->lpFirstFunction->value(alter) -
		this->lpSecondFunction->value(alter);
}

}
