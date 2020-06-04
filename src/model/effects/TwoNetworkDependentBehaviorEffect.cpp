/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TwoNetworkDependentBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * TwoNetworkDependentBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>
//#include "R_ext/Print.h"
#include "TwoNetworkDependentBehaviorEffect.h"
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
TwoNetworkDependentBehaviorEffect::TwoNetworkDependentBehaviorEffect(
	const EffectInfo * pEffectInfo) : BehaviorEffect(pEffectInfo)
{
	this->lfirstTotalAlterValues = 0;
	this->lfirstTotalInAlterValues = 0;
}

/**
 * Deallocates this effect object;
 */
TwoNetworkDependentBehaviorEffect::~TwoNetworkDependentBehaviorEffect()
{
	delete [] this->lfirstTotalAlterValues;
	delete [] this->lfirstTotalInAlterValues;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void TwoNetworkDependentBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	BehaviorEffect::initialize(pData, pState, period, pCache);
	string networkName1 = this->pEffectInfo()->interactionName1();
	string networkName2 = this->pEffectInfo()->interactionName2();

	this->lpFirstNetwork = pState->pNetwork(networkName1);
	this->lpSecondNetwork = pState->pNetwork(networkName2);

	if (!this->lpFirstNetwork)
	{
		throw logic_error("Network '" + networkName1 + "' expected.");
	}
	if (!this->lpSecondNetwork)
	{
		throw logic_error("Network '" + networkName2 + "' expected.");
	}
	if (this->lfirstTotalAlterValues)
	{
		delete [] this->lfirstTotalAlterValues;
	}
	if (this->lfirstTotalInAlterValues)
	{
		delete [] this->lfirstTotalInAlterValues;
	}
	this->lfirstTotalAlterValues = new double[this->lpFirstNetwork->n()];
	this->lfirstTotalInAlterValues = new double[this->lpFirstNetwork->m()];
}

/**
 * Returns the total alter covariate value for the given actor.
 */
double TwoNetworkDependentBehaviorEffect::firstTotalAlterValue(int i) const
{
	return this->lfirstTotalAlterValues[i];
}

/**
 * Returns the total in-alter covariate value for the given actor.
 */
double TwoNetworkDependentBehaviorEffect::firstTotalInAlterValue(int i) const
{
	return this->lfirstTotalInAlterValues[i];
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void TwoNetworkDependentBehaviorEffect::preprocessEgo(int ego)
{
	// set up the covariate based on current values of the network and behavior
	const Network * pFirstNetwork = this->pFirstNetwork();

	for (int i = 0; i < pFirstNetwork->n(); i++)
	{
		this->lfirstTotalAlterValues[i] = 0;
		if (pFirstNetwork->outDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pFirstNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->lfirstTotalAlterValues[i] += this->centeredValue(j);
// 				Rprintf("%d %f %d %d %d %d\n",
// 					j,
// 					this->centeredValue(j),
// 					this->period(),
			}
		}
		else
		{
			this->lfirstTotalAlterValues[i] = 0;
		}
//		Rprintf("%d %f\n", i,this->ltotalAlterValues[i]);
	}

	for (int i = 0; i < pFirstNetwork->m(); i++)
	{
		this->lfirstTotalInAlterValues[i] = 0;
		if (pFirstNetwork->inDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pFirstNetwork->inTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->lfirstTotalInAlterValues[i] += this->centeredValue(j);
			}
		}
		else
		{
			this->lfirstTotalInAlterValues[i] = 0;
		}
	}
}

}
