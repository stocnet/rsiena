/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConfigurationTable.cpp
 *
 * Description: This file contains the implementation of the ConfigurationTable
 * class.
 *****************************************************************************/

#include "ConfigurationTable.h"
#include "network/Network.h"
#include "model/tables/NetworkCache.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors, destructor, and initializers
// ----------------------------------------------------------------------------

/**
 * Creates an empty configuration table for the given network cache object.
 */
ConfigurationTable::ConfigurationTable(NetworkCache * pOwner)
{
	this->lpOwner = pOwner;
	this->lpNetwork = pOwner->pNetwork();
	this->ltable = new int[this->lpNetwork->n()];

	// This will make sure that the table is calculated on the first
	// call to the get(...) method, as the network modification count
	// starts with 0.

	this->llastModificationCount = -1;
}


/**
 * Deallocates this configuration table.
 */
ConfigurationTable::~ConfigurationTable()
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
int ConfigurationTable::get(int i)
{
	// If the network has changed, recalculate the table

	if (this->lpNetwork->modificationCount() != this->llastModificationCount)
	{
		this->calculate();
		this->llastModificationCount = this->lpNetwork->modificationCount();
	}

	return this->ltable[i];
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the network cache object owning this configuration table.
 */
NetworkCache * ConfigurationTable::pOwner() const
{
	return this->lpOwner;
}


/**
 * Returns the network this configuration table is associated with.
 */
const Network * ConfigurationTable::pNetwork() const
{
	return this->lpNetwork;
}


// ----------------------------------------------------------------------------
// Section: Protected methods
// ----------------------------------------------------------------------------

/**
 * Resets the internal array to zeroes.
 */
void ConfigurationTable::reset()
{
	for (int i = 0; i < this->lpNetwork->n(); i++)
	{
		this->ltable[i] = 0;
	}
}

}
