/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GenericNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GenericNetworkEffect.
 *****************************************************************************/

#include <stdexcept>
#include "GenericNetworkEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/effects/generic/AlterFunction.h"
#include "model/tables/Cache.h"
#include "model/EffectInfo.h"
#include "data/Data.h"

namespace siena
{

/**
 * Creates a new effect with the same function calculating the effect
 * and the statistic.
 */
GenericNetworkEffect::GenericNetworkEffect(const EffectInfo * pEffectInfo,
	AlterFunction * pFunction) : NetworkEffect(pEffectInfo)
{
	this->lpEffectFunction = pFunction;
	this->lpStatisticFunction = pFunction;
	this->lEffectType = pEffectInfo->effectType();
}


/**
 * Creates a new effect with different functions calculating the effect
 * and the statistic.
 */
GenericNetworkEffect::GenericNetworkEffect(const EffectInfo * pEffectInfo,
	AlterFunction * pEffectFunction,
	AlterFunction * pStatisticFunction) : NetworkEffect(pEffectInfo)
{
	this->lpEffectFunction = pEffectFunction;
	this->lpStatisticFunction = pStatisticFunction;
	this->lEffectType = pEffectInfo->effectType();
}


/**
 * Deallocates this effect.
 */
GenericNetworkEffect::~GenericNetworkEffect()
{
	if (this->lpEffectFunction == this->lpStatisticFunction)
	{
		delete this->lpEffectFunction;
	}
	else
	{
		delete this->lpEffectFunction;
		delete this->lpStatisticFunction;
	}
	this->lpEffectFunction = 0;
	this->lpStatisticFunction = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void GenericNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	this->lpEffectFunction->initialize(pData, pState, period, pCache);
	if (this->lpStatisticFunction != this->lpEffectFunction)
	{
		this->lpStatisticFunction->initialize(pData, pState, period, pCache);
	}
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] pSimulatedState the current simulated state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void GenericNetworkEffect::initialize(const Data * pData,
	State * pState, State * pSimulatedState,  
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, pSimulatedState, period, pCache);
	this->lpEffectFunction->initialize(pData, pState, pSimulatedState, period, pCache);
	if (this->lpStatisticFunction != this->lpEffectFunction)
	{
		this->lpStatisticFunction->initialize(pData, pState, pSimulatedState, period, pCache);
	}
// For the GMoM, only the estimation statistic is different.
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateContribution(...).
 */
void GenericNetworkEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);
	this->lpEffectFunction->preprocessEgo(ego);

	if (this->lpStatisticFunction != this->lpEffectFunction)
	{
		this->lpStatisticFunction->preprocessEgo(ego);
	}
}


/**
 * Assuming that the ego would flip the tie to the given actor,
 * this method calculates the change in the statistic corresponding
 * to this effect.
 */
double GenericNetworkEffect::calculateContribution(int alter) const
{
	return this->lpEffectFunction->value(alter);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double GenericNetworkEffect::tieStatistic(int alter)
{
	return this->lpStatisticFunction->value(alter);
}

}
