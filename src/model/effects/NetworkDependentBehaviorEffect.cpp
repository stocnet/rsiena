/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkDependentBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>
//#include "R_ext/Print.h"
#include "NetworkDependentBehaviorEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "network/IncidentTieIterator.h"
#include "model/State.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 */
NetworkDependentBehaviorEffect::NetworkDependentBehaviorEffect(
	const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo), //
	lSimulatedOffset(0), //
	ltotalAlterValues(0), //
	ltotalInAlterValues(0)
{
}

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 * @param simulatedState If `true` the value(), missing() and similarity()
 *        functions uses the simulated state, if any or the value at the end
 *        of the period.
 */
NetworkDependentBehaviorEffect::NetworkDependentBehaviorEffect(
	const EffectInfo * pEffectInfo, bool simulatedState) :
	BehaviorEffect(pEffectInfo), //
	lSimulatedOffset(simulatedState ? 1 : 0), //
	ltotalAlterValues(0), //
	ltotalInAlterValues(0)
{
}

/**
 * Deallocates this effect object;
 */
NetworkDependentBehaviorEffect::~NetworkDependentBehaviorEffect()
{
	delete [] this->ltotalAlterValues;
	delete [] this->ltotalInAlterValues;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void NetworkDependentBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	BehaviorEffect::initialize(pData, pState, period, pCache);
	string networkName = this->pEffectInfo()->interactionName1();
	this->lpNetwork = pState->pNetwork(networkName);

	if (!this->lpNetwork) {
		throw logic_error("Network '" + networkName + "' expected.");
	}

	// clear old value arrays
	if (this->ltotalAlterValues) delete [] this->ltotalAlterValues;
	if (this->ltotalInAlterValues) delete [] this->ltotalInAlterValues;
	// create new value arrays
	this->ltotalAlterValues = new double[this->lpNetwork->n()];
	this->ltotalInAlterValues = new double[this->lpNetwork->m()];
}

void NetworkDependentBehaviorEffect::initialize(const Data *pData,
	State *pState, State *pSimulatedState, int period, Cache *pCache)
{
	BehaviorEffect::initialize(pData, pState, period, pCache);
	string networkName = this->pEffectInfo()->interactionName1();
	this->lpNetwork = pState->pNetwork(networkName);

	if (!this->lpNetwork)
	{
		throw logic_error("Network '" + networkName + "' expected.");
	}

	// Select network state.
	if (lSimulatedOffset == 1)
	{
		this->lpNetwork = pSimulatedState->pNetwork(networkName);
	}
	else
	{
		this->lpNetwork = pState->pNetwork(networkName);
	}

	// clear old value arrays
	if (this->ltotalAlterValues) delete [] this->ltotalAlterValues;
	if (this->ltotalInAlterValues) delete [] this->ltotalInAlterValues;
	// create new value arrays
	this->ltotalAlterValues = new double[this->lpNetwork->n()];
	this->ltotalInAlterValues = new double[this->lpNetwork->m()];
}

/**
 * Returns the total alter covariate value for the given actor.
 */
double NetworkDependentBehaviorEffect::totalAlterValue(int i) const
{
	return this->ltotalAlterValues[i];
}

/**
 * Returns the total in-alter covariate value for the given actor.
 */
double NetworkDependentBehaviorEffect::totalInAlterValue(int i) const
{
	return this->ltotalInAlterValues[i];
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void NetworkDependentBehaviorEffect::preprocessEgo(int ego)
{
	// set up the covariate based on current values of the network and behavior
	const Network * pNetwork = this->pNetwork();

	for (int i = 0; i < pNetwork->n(); i++)
	{
		this->ltotalAlterValues[i] = 0;
		if (pNetwork->outDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->ltotalAlterValues[i] += this->centeredValue(j);
// 				Rprintf("%d %f %d %d %d %d\n",
// 					j, this->centeredValue(j), this->period(),
			}
		}
		else
		{
			this->ltotalAlterValues[i] = 0;
		}
//		Rprintf("%d %f\n", i,this->ltotalAlterValues[i]);
	}

	for (int i = 0; i < pNetwork->m(); i++)
	{
		this->ltotalInAlterValues[i] = 0;
		if (pNetwork->inDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->inTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->ltotalInAlterValues[i] += this->centeredValue(j);
			}
		}
		else
		{
			this->ltotalInAlterValues[i] = 0;
		}
	}
}

}
