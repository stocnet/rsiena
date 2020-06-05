/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkCache.cpp
 *
 * Description: This file contains the implementation of the class
 * NetworkCache.
 *****************************************************************************/

#include "NetworkCache.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "model/tables/TwoPathTable.h"
#include "model/tables/CriticalInStarTable.h"
#include "model/tables/BetweennessTable.h"

namespace siena
{

/**
 * Constructs a network cache for the given network.
 */
NetworkCache::NetworkCache(const Network * pNetwork)
{
	this->lpNetwork = pNetwork;

	this->loutTieValues = new int[pNetwork->m()];

	this->loneModeNetwork =
		dynamic_cast<const OneModeNetwork *>(pNetwork) != 0;

	// Many of the configurations can be seen
	// as special cases of generalized two-paths, where one can specify
	// the directions of the first and the second step. For instance, the
	// number of in-stars between i and j equals the number of possible ways
	// of reaching j, if we have to traverse one outgoing tie of actor i,
	// say (i,h), followed by on incoming tie of actor h.

	if (this->loneModeNetwork)
	{
		this->linTieValues = new int[pNetwork->n()];
		this->lpTwoPathTable =
			new TwoPathTable(this, FORWARD, FORWARD);
		this->lpReverseTwoPathTable =
			new TwoPathTable(this, BACKWARD, BACKWARD);
		this->lpOutStarTable =
			new TwoPathTable(this, BACKWARD, FORWARD);
		this->lpCriticalInStarTable = new CriticalInStarTable(this);
		this->lpRRTable =
			new TwoPathTable(this, RECIPROCAL, RECIPROCAL);
		this->lpRFTable =
			new TwoPathTable(this, RECIPROCAL, FORWARD);
		this->lpRBTable =
			new TwoPathTable(this, RECIPROCAL, BACKWARD);
		this->lpFRTable =
			new TwoPathTable(this, FORWARD, RECIPROCAL);
		this->lpBRTable =
			new TwoPathTable(this, BACKWARD, RECIPROCAL);
		this->lpBetweennessTable = new BetweennessTable(this);
	}
	else
	{
		this->linTieValues = 0;
		this->lpTwoPathTable = 0;
		this->lpReverseTwoPathTable = 0;
		this->lpOutStarTable = 0;
		this->lpCriticalInStarTable = 0;
		this->lpRRTable = 0;
		this->lpRFTable = 0;
		this->lpRBTable = 0;
		this->lpFRTable = 0;
		this->lpBRTable = 0;
		this->lpBetweennessTable = 0;
	}

	this->lpInStarTable = new TwoPathTable(this, FORWARD, BACKWARD);

	this->initialize(-1);
}


/**
 * Destroys this network cache object.
 */
NetworkCache::~NetworkCache()
{
	delete[] this->loutTieValues;
	delete[] this->linTieValues;
	delete this->lpTwoPathTable;
	delete this->lpReverseTwoPathTable;
	delete this->lpOutStarTable;
	delete this->lpCriticalInStarTable;
	delete this->lpRRTable;
	delete this->lpRFTable;
	delete this->lpRBTable;
	delete this->lpFRTable;
	delete this->lpBRTable;
	delete this->lpBetweennessTable;
	delete this->lpInStarTable;

	this->loutTieValues = 0;
	this->linTieValues = 0;
	this->lpTwoPathTable = 0;
	this->lpReverseTwoPathTable = 0;
	this->lpOutStarTable = 0;
	this->lpCriticalInStarTable = 0;
	this->lpRRTable = 0;
	this->lpRFTable = 0;
	this->lpRBTable = 0;
	this->lpFRTable = 0;
	this->lpBRTable = 0;
	this->lpBetweennessTable = 0;
	this->lpInStarTable = 0;
}


void NetworkCache::initialize(int ego)
{
	// Out-tie indicators

	for (int i = 0; i < this->lpNetwork->m(); i++)
	{
		this->loutTieValues[i] = 0;
	}

	if (ego >= 0 && ego < this->lpNetwork->n())
	{
		for (IncidentTieIterator iter = this->lpNetwork->outTies(ego);
			iter.valid();
			iter.next())
		{
			this->loutTieValues[iter.actor()] = iter.value();
		}
	}

	// In-tie indicators

	if (this->loneModeNetwork)
	{
		for (int i = 0; i < this->lpNetwork->n(); i++)
		{
			this->linTieValues[i] = 0;
		}

		if (ego >= 0 && ego < this->lpNetwork->n())
		{
			for (IncidentTieIterator iter = this->lpNetwork->inTies(ego);
				iter.valid();
				iter.next())
			{
				this->linTieValues[iter.actor()] = iter.value();
			}
		}
	}

	// Initialize all configuration tables

	if (this->loneModeNetwork)
	{
		this->lpTwoPathTable->initialize(ego);
		this->lpReverseTwoPathTable->initialize(ego);
		this->lpOutStarTable->initialize(ego);
		this->lpCriticalInStarTable->initialize(ego);
		this->lpRRTable->initialize(ego);
		this->lpRFTable->initialize(ego);
		this->lpRBTable->initialize(ego);
		this->lpFRTable->initialize(ego);
		this->lpBRTable->initialize(ego);
	}

	this->lpInStarTable->initialize(ego);
}

}
