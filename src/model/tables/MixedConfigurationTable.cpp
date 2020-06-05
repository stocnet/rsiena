/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedConfigurationTable.cpp
 *
 * Description: This file contains the implementation of the
 * MixedConfigurationTable class.
 *****************************************************************************/

#include "MixedConfigurationTable.h"
#include "network/Network.h"
#include "model/tables/TwoNetworkCache.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors, destructor, and initializers
// ----------------------------------------------------------------------------

/**
 * Creates an empty configuration table for the given two network cache object.
 */
MixedConfigurationTable::MixedConfigurationTable(TwoNetworkCache * pOwner)
{
	this->lpOwner = pOwner;
	this->lpFirstNetwork = pOwner->pFirstNetwork();
	this->lpSecondNetwork = pOwner->pSecondNetwork();

	// Make sure the table is large enough for any of the involved sets
	// of actors.

	this->ltableSize =
		std::max(std::max(this->lpFirstNetwork->n(), this->lpFirstNetwork->m()),
			std::max(this->lpSecondNetwork->n(), this->lpSecondNetwork->m()));

	this->ltable = new int[this->ltableSize];

	// This will make sure that the table is calculated on the first
	// call to the get(...) method, as the network modification count
	// starts with 0.

	this->llastFirstModificationCount = -1;
	this->llastSecondModificationCount = -1;
}


/**
 * Deallocates this configuration table.
 */
MixedConfigurationTable::~MixedConfigurationTable()
{
	delete[] this->ltable;
	this->ltable = 0;
}


// ----------------------------------------------------------------------------
// Section: Public interface
// ----------------------------------------------------------------------------

/**
 * Returns the number of configurations corresponding to the given actor.
 */
int MixedConfigurationTable::get(int i)
{
	// If the network has changed, recalculate the table

	if ((this->lpFirstNetwork->modificationCount() !=
		this->llastFirstModificationCount) ||
		(this->lpSecondNetwork->modificationCount() !=
			this->llastSecondModificationCount))
	{
		this->calculate();
		this->llastFirstModificationCount =
			this->lpFirstNetwork->modificationCount();
		this->llastSecondModificationCount =
			this->lpSecondNetwork->modificationCount();
	}

	return this->ltable[i];
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the network cache object owning this configuration table.
 */
TwoNetworkCache * MixedConfigurationTable::pOwner() const
{
	return this->lpOwner;
}


/**
 * Returns the first network this configuration table is associated with.
 */
const Network * MixedConfigurationTable::pFirstNetwork() const
{
	return this->lpFirstNetwork;
}


/**
 * Returns the second network this configuration table is associated with.
 */
const Network * MixedConfigurationTable::pSecondNetwork() const
{
	return this->lpSecondNetwork;
}

// ----------------------------------------------------------------------------
// Section: Protected methods
// ----------------------------------------------------------------------------

/**
 * Resets the internal array to zeroes.
 */
void MixedConfigurationTable::reset()
{
	for (int i = 0; i < this->ltableSize; i++)
	{
		this->ltable[i] = 0;
	}
}

}
