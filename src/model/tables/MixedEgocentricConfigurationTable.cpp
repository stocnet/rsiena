/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedEgocentricConfigurationTable.h
 *
 * Description: This file defines the class MixedEgocentriConfigurationTable.
 *****************************************************************************/

#include "MixedEgocentricConfigurationTable.h"

namespace siena
{

MixedEgocentricConfigurationTable::MixedEgocentricConfigurationTable(
	TwoNetworkCache * pOwner) : MixedConfigurationTable(pOwner)
{
	// The table has not been calculated for any ego yet.

	this->lego = -1;
	this->lupdated = false;
}


MixedEgocentricConfigurationTable::~MixedEgocentricConfigurationTable()
{
}


/**
 * Initializes the table for the given ego.
 */
void MixedEgocentricConfigurationTable::initialize(int ego)
{
	// Store the ego
	this->lego = ego;

	// Make sure the table is recalculated in the next call to get(...)
	this->lupdated = false;
}


/**
 * Returns the number of configurations corresponding to the given actor.
 */
int MixedEgocentricConfigurationTable::get(int i)
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
int MixedEgocentricConfigurationTable::ego() const
{
	return this->lego;
}

}
