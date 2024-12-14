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
	ltotalInAlterValues(0),
	lnumberAlterHigher(0),
	lnumberAlterLower(0),
	lnumberAlterEqual(0),
	lnumberAlterHigherPop(0),
	lnumberAlterLowerPop(0),
	lnumberAlterEqualPop(0)
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
	ltotalInAlterValues(0),
	lnumberAlterHigher(0),
	lnumberAlterLower(0),
	lnumberAlterEqual(0),
	lnumberAlterHigherPop(0),
	lnumberAlterLowerPop(0),
	lnumberAlterEqualPop(0)
{
}

/**
 * Deallocates this effect object;
 */
NetworkDependentBehaviorEffect::~NetworkDependentBehaviorEffect()
{
	delete [] this->ltotalAlterValues;
	delete [] this->ltotalInAlterValues;
	delete [] this->lnumberAlterHigher;
	delete [] this->lnumberAlterLower;
	delete [] this->lnumberAlterEqual;
	delete [] this->lnumberAlterHigherPop;
	delete [] this->lnumberAlterLowerPop;
	delete [] this->lnumberAlterEqualPop;
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
	if (this->ltotalAlterValues)  delete [] this->ltotalAlterValues;
	if (this->ltotalInAlterValues) delete [] this->ltotalInAlterValues;
	if (this->lnumberAlterHigher) delete [] this->lnumberAlterHigher;
	if (this->lnumberAlterLower)  delete [] this->lnumberAlterLower;
	if (this->lnumberAlterEqual)  delete [] this->lnumberAlterEqual;
	if (this->lnumberAlterHigherPop) delete [] this->lnumberAlterHigherPop;
	if (this->lnumberAlterLowerPop)  delete [] this->lnumberAlterLowerPop;
	if (this->lnumberAlterEqualPop)  delete [] this->lnumberAlterEqualPop;
	// create new value arrays
	this->ltotalAlterValues  = new double[this->lpNetwork->n()];
	this->ltotalInAlterValues = new double[this->lpNetwork->m()];
	this->lnumberAlterHigher = new int[this->lpNetwork->n()];
	this->lnumberAlterLower  = new int[this->lpNetwork->n()];
	this->lnumberAlterEqual  = new int[this->lpNetwork->n()];
	this->lnumberAlterHigherPop = new int[this->lpNetwork->n()];
	this->lnumberAlterLowerPop  = new int[this->lpNetwork->n()];
	this->lnumberAlterEqualPop  = new int[this->lpNetwork->n()];
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
	if (this->ltotalAlterValues)  delete [] this->ltotalAlterValues;
	if (this->ltotalInAlterValues) delete [] this->ltotalInAlterValues;
	if (this->lnumberAlterEqual)  delete [] this->lnumberAlterEqual;
	if (this->lnumberAlterHigher) delete [] this->lnumberAlterHigher;
	if (this->lnumberAlterLower)  delete [] this->lnumberAlterLower;
	if (this->lnumberAlterEqualPop)  delete [] this->lnumberAlterEqualPop;
	if (this->lnumberAlterHigherPop) delete [] this->lnumberAlterHigherPop;
	if (this->lnumberAlterLowerPop)  delete [] this->lnumberAlterLowerPop;
	// create new value arrays
	this->ltotalAlterValues  = new double[this->lpNetwork->n()];
	this->ltotalInAlterValues = new double[this->lpNetwork->m()];
	this->lnumberAlterEqual  = new int[this->lpNetwork->n()];
	this->lnumberAlterHigher = new int[this->lpNetwork->n()];
	this->lnumberAlterLower  = new int[this->lpNetwork->n()];
	this->lnumberAlterEqualPop  = new int[this->lpNetwork->n()];
	this->lnumberAlterHigherPop = new int[this->lpNetwork->n()];
	this->lnumberAlterLowerPop  = new int[this->lpNetwork->n()];
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
 * Returns the number of higher alter covariate values for the given actor.
 */
int NetworkDependentBehaviorEffect::numberAlterHigher(int i) const
{
	return this->lnumberAlterHigher[i];
}

/**
 * Returns the number of lower alter covariate values for the given actor.
 */
int NetworkDependentBehaviorEffect::numberAlterLower(int i) const
{
	return this->lnumberAlterLower[i];
}


/**
 * Returns the number of equal alter covariate values for the given actor.
 */
int NetworkDependentBehaviorEffect::numberAlterEqual(int i) const
{
	return this->lnumberAlterEqual[i];
}



/**
 * Returns the number of higher alter covariate values,
 *  weighted by alter indegree, for the given actor.
 */
int NetworkDependentBehaviorEffect::numberAlterHigherPop(int i) const
{
	return this->lnumberAlterHigherPop[i];
}

/**
 * Returns the number of lower alter covariate values,
 *  weighted by alter indegree, for the given actor.
 */
int NetworkDependentBehaviorEffect::numberAlterLowerPop(int i) const
{
	return this->lnumberAlterLowerPop[i];
}


/**
 * Returns the number of equal alter covariate values,
 *  weighted by alter indegree, for the given actor.
 */
int NetworkDependentBehaviorEffect::numberAlterEqualPop(int i) const
{
	return this->lnumberAlterEqualPop[i];
}

/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void NetworkDependentBehaviorEffect::preprocessEgo(int ego)
{
	BehaviorEffect::preprocessEgo(ego);
	// set up the covariate based on current values of the network and behavior
	const Network * pNetwork = this->pNetwork();

	for (int i = 0; i < pNetwork->n(); i++)
	{
		int vego = this->value(i); // non-centered
		this->ltotalAlterValues[i] = 0;
		this->lnumberAlterHigher[i] = 0;
		this->lnumberAlterLower[i] = 0;
		this->lnumberAlterEqual[i] = 0;
		this->lnumberAlterHigherPop[i] = 0;
		this->lnumberAlterLowerPop[i] = 0;
		this->lnumberAlterEqualPop[i] = 0;
		if (pNetwork->outDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->ltotalAlterValues[i] += this->centeredValue(j);
				if (this->value(j) > vego)
				{
					lnumberAlterHigher[i]++;
					lnumberAlterHigherPop[i] += pNetwork->inDegree(j);
				}
				else
				{
					if (this->value(j) < vego)
					{
						lnumberAlterLower[i]++;
						lnumberAlterLowerPop[i] += pNetwork->inDegree(j);
					}
					else
					{
						lnumberAlterEqual[i]++;
						lnumberAlterHigherPop[i] += pNetwork->inDegree(j);
					}
				}
			}
// 				Rprintf("%d %f %d %d %d %d\n",
// 					j, this->centeredValue(j), this->period(),
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
	}
}

}
