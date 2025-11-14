/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: https://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EqualsZeroFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * EqualsZeroFunction.
 * Testing of equality with a different constant value 
 * could be easily accomodated.
 *****************************************************************************/

#include <cmath>
#include "EqualsZeroFunction.h"

namespace siena
{

EqualsZeroFunction::EqualsZeroFunction(AlterFunction * pFunction, bool minus)
{
	this->lpFunction = pFunction;
	this->lminus = minus;
}


/**
 * Deallocates this function.
 */
EqualsZeroFunction::~EqualsZeroFunction()
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
void EqualsZeroFunction::initialize(const Data * pData,
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
 * calling EqualsZeroFunction::value(...).
 */
void EqualsZeroFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	this->lpFunction->preprocessEgo(ego);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double EqualsZeroFunction::value(int alter) const
{
	double statistic = 0;
	if (fabs(this->lpFunction->value(alter)) < EPSILON)
	{		
		if (this->lminus)
		{
			statistic = -1;			
		}
		else
		{
			statistic = 1;			
		}
	}
	return statistic;
}

}
