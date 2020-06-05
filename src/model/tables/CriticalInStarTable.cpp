/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CriticalInStarTable.cpp
 *
 * Description: This file contains the implementation of the class
 * CriticalInStarTable.
 *****************************************************************************/

#include "CriticalInStarTable.h"
#include "model/tables/NetworkCache.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"
#include "network/Network.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Creates a new table for storing critical in-stars.
 */
CriticalInStarTable::CriticalInStarTable(NetworkCache * pOwner) :
	EgocentricConfigurationTable(pOwner)
{
}


// ----------------------------------------------------------------------------
// Section: ConfigurationTable implementation
// ----------------------------------------------------------------------------

/**
 * Calculates the number of critical in-stars between the ego and all
 * other actors.
 */
void CriticalInStarTable::calculate()
{
	// Reset the counters to zeroes
	this->reset();

	// We will need the two-path counts
	ConfigurationTable * pTwoPathTable = this->pOwner()->pTwoPathTable();

	const Network * pNetwork = this->pNetwork();

	// Consider each outgoing tie of the ego in turn.

	for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
		iter.valid();
		iter.next())
	{
		// Get the receiver of the outgoing tie.
		int h = iter.actor();

		// If there are more than one two-paths between i and h, then no
		// in-star pointing to h can be critical.

		if (pTwoPathTable->get(h) == 0)
		{
			// If there are no two-paths between i and h, then every in-star
			// pointing to h is critical. Just iterate over incoming ties of
			// h to reach them all.

			for (IncidentTieIterator iter1 = pNetwork->inTies(h);
				iter1.valid();
				iter1.next())
			{
				this->ltable[iter1.actor()]++;
			}
		}
		else if (pTwoPathTable->get(h) == 1)
		{
			// If there is exactly one two-path between i and h, say
			// i -> j -> h, then <(i,h), (j,h)> is the only critical in-star
			// pointing to h. Just iterate over incoming ties of h to
			// find the actor j.

			// This can be optimized by storing one of the intermediary
			// actors between ego and each of the other actors when
			// counting two-paths. The the only intermediary actor could
			// be accessed in constant time here.

			bool found = false;

			for (IncidentTieIterator iter1 = pNetwork->inTies(h);
				iter1.valid() && !found;
				iter1.next())
			{
				int j = iter1.actor();

				if (this->pOwner()->outTieExists(j))
				{
					this->ltable[j]++;
					found = true;
				}
			}
		}
	}
}

}
