/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkDependentContinuousEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "NetworkDependentContinuousEffect.h"
#include "model/State.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 */
NetworkDependentContinuousEffect::NetworkDependentContinuousEffect(
	const EffectInfo * pEffectInfo) : ContinuousEffect(pEffectInfo)
{
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void NetworkDependentContinuousEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	ContinuousEffect::initialize(pData, pState, period, pCache);
	string networkName = this->pEffectInfo()->interactionName1();

	this->lpNetwork = pState->pNetwork(networkName);

	if (!this->lpNetwork)
	{
		throw logic_error("Network '" + networkName + "' expected.");
	}
}

}
