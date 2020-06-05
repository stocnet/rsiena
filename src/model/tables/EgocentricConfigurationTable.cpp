#include "EgocentricConfigurationTable.h"

namespace siena
{

EgocentricConfigurationTable::EgocentricConfigurationTable(
	NetworkCache * pOwner) : ConfigurationTable(pOwner)
{
	// The table has not been calculated for any ego yet.

	this->lego = -1;
	this->lupdated = false;
}


EgocentricConfigurationTable::~EgocentricConfigurationTable()
{
}


/**
 * Initializes the table for the given ego.
 */
void EgocentricConfigurationTable::initialize(int ego)
{
	// Store the ego
	this->lego = ego;

	// Make sure the table is recalculated in the next call to get(...)
	this->lupdated = false;
}


/**
 * Returns the number of configurations corresponding to the given actor.
 */
int EgocentricConfigurationTable::get(int i)
{
	if (!this->lupdated)
	{
		this->calculate();
		this->lupdated = true;
	}

	return this->ltable[i];
}


/**
 * Returns the current ego of this configuration table.
 */
int EgocentricConfigurationTable::ego() const
{
	return this->lego;
}

}
