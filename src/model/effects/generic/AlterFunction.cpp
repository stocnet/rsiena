/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * AlterFunction.
 *****************************************************************************/

#include "AlterFunction.h"
#include "data/Data.h"
#include "model/State.h"
#include "model/tables/Cache.h"

namespace siena
{

AlterFunction::AlterFunction()
{
	this->lego = -1;
}

AlterFunction::~AlterFunction()
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void AlterFunction::initialize(const Data * pData,
		State * pState, int period, Cache * pCache)
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] pSimulatedState the current simulated state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void AlterFunction::initialize(const Data * pData,
	State * pState, State * pSimulatedState,  
	int period,
	Cache * pCache)
{
}

/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling AlterFunction::value(...).
 */
void AlterFunction::preprocessEgo(int ego)
{
	this->lego = ego;
}

}
