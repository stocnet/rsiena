/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BetweennessTable.cpp
 *
 * Description: This file contains the implementation of the BetweennessTable
 * class.
 *****************************************************************************/

#include "BetweennessTable.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Constructor.
 */
BetweennessTable::BetweennessTable(NetworkCache * pOwner) :
	ConfigurationTable(pOwner)
{
}


// ----------------------------------------------------------------------------
// Section: ConfigurationTable implementation
// ----------------------------------------------------------------------------

/**
 * Calculates the betweenness counts for all actors.
 */
void BetweennessTable::calculate()
{
	// Reset the counters to zeroes
	this->reset();

	// One-mode network is assumed

	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();

	// Helper array
	int * mark = new int[n];

	for (int i = 0; i < n; i++)
	{
		mark[i] = -1;
	}

	for (int i = 0; i < n; i++)
	{
		// Mark the out-neighbors of actor i

		for (IncidentTieIterator iter = pNetwork->outTies(i);
			iter.valid();
			iter.next())
		{
			mark[iter.actor()] = i;
		}

		// Consider each two-path starting at i

		for (IncidentTieIterator iterI = pNetwork->outTies(i);
			iterI.valid();
			iterI.next())
		{
			int j = iterI.actor();

			for (IncidentTieIterator iterJ = pNetwork->outTies(j);
				iterJ.valid();
				iterJ.next())
			{
				int h = iterJ.actor();

				if (i != h && mark[h] != i)
				{
					// We have found a two-path i -> j -> h with no tie
					// from i to h, so we increase the betweenness of j.

					this->ltable[j]++;
				}
			}
		}
	}

	delete[] mark;
}

}
