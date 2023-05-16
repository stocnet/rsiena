/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkWithPrimaryEffect.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkWithPrimaryEffect class.
 *****************************************************************************/

#include "NetworkWithPrimaryEffect.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"
#include "network/IncidentTieIterator.h"
#include "network/UnionNeighborIterator.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
NetworkWithPrimaryEffect::NetworkWithPrimaryEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
	this->lprimary = 0;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void NetworkWithPrimaryEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	int n = this->pNetwork()->n();
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());
	if (!pONetwork)
	{
		throw runtime_error(
			"One-mode network expected in NetworkWithPrimaryEffect");
	}

	delete [] this->lprimary;
	this->lprimary = new bool[n];
}

/**
 * Destructor.
 */
NetworkWithPrimaryEffect::~NetworkWithPrimaryEffect()
{
	delete [] this->lprimary;
}


/**
 * Calculates primary setting properties for ego:
 * lprimary, lprimDegree, llogNonPrimary
 * also regular outdegree ld
 */
void NetworkWithPrimaryEffect::primaryProperties(int ego, const Network * pNetwork)
{
	int n = pNetwork->n();
	const OneModeNetwork * pOneModeNetwork =  
		dynamic_cast<const OneModeNetwork *>(pNetwork);
	this->lprimDegree = 0;	
	for (int i = 0; i < n; i++)
	{
		this->lprimary[i] = false;
	}
	
	for (UnionNeighborIterator iter = pOneModeNetwork->eitherTies(ego);
		iter.valid();
		iter.next())
	{
		int j = iter.actor();
		if (!this->lprimary[j])			
		{
			this->lprimDegree++;
			lprimary[j] = true;
		}
		for (UnionNeighborIterator iter2 = pOneModeNetwork->eitherTies(j);
			iter2.valid();
			iter2.next())
		{
			int k = iter2.actor();
			if (!this->lprimary[k])			
			{
				this->lprimDegree++;
				this->lprimary[k] = true;
			}
		}
	}
	if (this->lprimary[ego])
	{		
		this->lprimDegree--;
		this->lprimary[ego] = false;
	}
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void NetworkWithPrimaryEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);
	const Network * pNetwork = this->pNetwork();
	this->primaryProperties(ego, pNetwork);
}


int  NetworkWithPrimaryEffect::primaryDegree() const
/**
 * degree of ego in primary setting, after preprocessEgo(int ego) 
 */
{
	return this->lprimDegree;	 
}

bool NetworkWithPrimaryEffect::inPrimarySetting(int alter) const
/**
 * whether alter is in ego's in primary setting, after preprocessEgo(int ego) 
 */
{
	return this->lprimary[alter];	 
}

}
