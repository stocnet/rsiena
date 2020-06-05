/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TwoNetworkCache.cpp
 *
 * Description: This file contains the implementation of the class
 * TwoNetworkCache.
 *****************************************************************************/

#include "TwoNetworkCache.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "model/tables/MixedTwoPathTable.h"

namespace siena
{

/**
 * Constructs a  twonetwork cache for the given networks.
 * Note that first network = W (explanatory network),
 * second network = X (dependent network)
 */
TwoNetworkCache::TwoNetworkCache(const Network * pFirstNetwork,
	const Network * pSecondNetwork)
{
	this->lpFirstNetwork = pFirstNetwork;
	this->lpSecondNetwork = pSecondNetwork;

	this->lfirstOutTieValues = new int[pFirstNetwork->m()];
	this->lsecondOutTieValues = new int[pSecondNetwork->m()];

	this->loneModeFirstNetwork =
		dynamic_cast<const OneModeNetwork *>(pFirstNetwork) != 0;
	this->loneModeSecondNetwork =
		dynamic_cast<const OneModeNetwork *>(pSecondNetwork) != 0;

	// Many of the configurations can be seen
	// as special cases of generalized two-paths, where one can specify
	// the directions of the first and the second step. For instance, the
	// number of in-stars between i and j equals the number of possible ways
	// of reaching j, if we have to traverse one outgoing tie of actor i,
	// say (i,h), followed by on incoming tie of actor h.

	// Not sure how many are relevant to W (first) a one mode, X (second)
	// a one or twomode.

	if (this->loneModeFirstNetwork)
	{
		this->lpTwoPathTable =
			new MixedTwoPathTable(this, FORWARD, FORWARD);
	}
	else
	{
		this->lpTwoPathTable = 0;
	}

	this->lpInStarTable =
		new MixedTwoPathTable(this, FORWARD, BACKWARD);

	this->lpOutStarTable =
		new MixedTwoPathTable(this, BACKWARD, FORWARD);

	this->initialize(-1);
}


/**
 * Destroys this network cache object.
 */
TwoNetworkCache::~TwoNetworkCache()
{
	delete[] this->lfirstOutTieValues;
	delete[] this->lsecondOutTieValues;
	delete this->lpTwoPathTable;
	delete this->lpInStarTable;
	delete this->lpOutStarTable;

	this->lfirstOutTieValues = 0;
	this->lsecondOutTieValues = 0;
	this->lpTwoPathTable = 0;
	this->lpInStarTable = 0;
	this->lpOutStarTable = 0;
}


void TwoNetworkCache::initialize(int ego)
{
	// Out-tie indicators

	for (int i = 0; i < this->lpFirstNetwork->m(); i++)
	{
		this->lfirstOutTieValues[i] = 0;
	}

	if (ego >= 0 && ego < this->lpFirstNetwork->n())
	{
		for (IncidentTieIterator iter = this->lpFirstNetwork->outTies(ego);
			iter.valid();
			iter.next())
		{
			this->lfirstOutTieValues[iter.actor()] = iter.value();
		}
	}

	for (int i = 0; i < this->lpSecondNetwork->m(); i++)
	{
		this->lsecondOutTieValues[i] = 0;
	}

	if (ego >= 0 && ego < this->lpSecondNetwork->n())
	{
		this->lsecondOutDegree = 0;
		for (IncidentTieIterator iter = this->lpSecondNetwork->outTies(ego);
			iter.valid();
			iter.next())
		{
			this->lsecondOutTieValues[iter.actor()] = iter.value();
			this->lsecondOutDegree++;
		}
	}


	// Initialize all configuration tables

	if (this->loneModeFirstNetwork)
	{
		this->lpTwoPathTable->initialize(ego);
	}

	this->lpInStarTable->initialize(ego);

	this->lpOutStarTable->initialize(ego);

}

}
