/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Effect.cpp
 *
 * Description: This file contains the implementation of the class Effect.
 *****************************************************************************/

#include "Effect.h"
#include "model/EffectInfo.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction, destruction
// ----------------------------------------------------------------------------

/**
 * Constructor.
 */
Effect::Effect(const EffectInfo * pEffectInfo)
{
	this->lpEffectInfo = pEffectInfo;
	this->lparameter = pEffectInfo->parameter();
	this->lperiod = 0;
}


/**
 * Destroys this effect.
 */
Effect::~Effect()
{
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the effect info object this effect is based on.
 */
const EffectInfo * Effect::pEffectInfo() const
{
	return this->lpEffectInfo;
}


// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void Effect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	this->lperiod = period;
	this->lpCache = pCache;
}

void Effect::parameter(double value)
{
	this->lparameter = value;
}
}
