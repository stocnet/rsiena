/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: Function.cpp
 * 
 * Description: This file implements the Function class.
 *****************************************************************************/

#include "Function.h"
#include "utils/Utils.h"
#include "model/effects/Effect.h"


namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors, destructors.
// ----------------------------------------------------------------------------

/**
 * Deallocates this function.
 */
Function::~Function()
{
	// Delete the effects stored in this function.
	deallocateVector(this->leffects);
}


// ----------------------------------------------------------------------------
// Section: Effect management
// ----------------------------------------------------------------------------

/**
 * Adds the given effect to this function.
 */
void Function::addEffect(Effect * pEffect)
{
	this->leffects.push_back(pEffect);
}


/**
 * Returns the collection of effects defining this function.
 */
const std::vector<Effect *> & Function::rEffects() const
{
	return this->leffects;
}

}
