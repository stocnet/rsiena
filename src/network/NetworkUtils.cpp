/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DataUtils.cpp
 *
 * Description: This module contains some utilities specific to the
 * 'data' library.
 *****************************************************************************/
#include <vector>
#include <set>
#include "NetworkUtils.h"
#include "IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "network/TieIterator.h"
#include "network/Network.h"

namespace siena
{

/**
 * Returns the number of actors iterated over by both of the given iterators.
 */
int commonActorCount(IncidentTieIterator iter1, IncidentTieIterator iter2)
{
	// Note that ties incident to an actor are sorted in an increasing order
	// of its neighbors.

	int count = 0;

	CommonNeighborIterator iter(iter1, iter2);

	while (iter.valid())
	{
		count++;
		iter.next();
	}

	return count;
}


/**
 * Creates a new network representing the symmetric difference of the
 * two given networks.
 */
Network * symmetricDifference(const Network * pNetwork1,
	const Network * pNetwork2)
{
	Network * pDifference = new Network(pNetwork1->n(), pNetwork1->m());

	for (int i = 0; i < pNetwork1->n(); i++)
	{
		IncidentTieIterator iter1 = pNetwork1->outTies(i);
		IncidentTieIterator iter2 = pNetwork2->outTies(i);

		while (iter1.valid() && iter2.valid())
		{
			if (iter1.actor() < iter2.actor())
			{
				pDifference->setTieValue(i, iter1.actor(), 1);
				iter1.next();
			}
			else if (iter1.actor() > iter2.actor())
			{
				pDifference->setTieValue(i, iter2.actor(), 1);
				iter2.next();
			}
			else
			{
				iter1.next();
				iter2.next();
			}
		}

		while (iter1.valid())
		{
			pDifference->setTieValue(i, iter1.actor(), 1);
			iter1.next();
		}

		while (iter2.valid())
		{
			pDifference->setTieValue(i, iter2.actor(), 1);
			iter2.next();
		}
	}

	return pDifference;
}

/**
 * Creates a network by removing from the first network the ties which are
 * present in the second network.
 */
void subtractNetwork(Network * pNetwork,
	const Network * pSubtrahendNetwork)
{
	for (TieIterator iter = pSubtrahendNetwork->ties();
		iter.valid();
		iter.next())
	{
		pNetwork->setTieValue(iter.ego(), iter.alter(), 0);
	}
}

/**
 * Replaces the values of the first network with values in the second network
 * for ties that are present in the third network.
 */
void replaceNetwork(Network * pNetwork,
	const Network * pValueNetwork, const Network * pDecisionNetwork )
{
	for (TieIterator iter = pDecisionNetwork->ties();
		iter.valid();
		iter.next())
	{
		pNetwork->setTieValue(iter.ego(), iter.alter(),
			pValueNetwork->tieValue(iter.ego(), iter.alter()));
	}
}

/**
 * Create primary setting based on a network
 */
std::vector<int> * primarySetting(const Network * pNetwork, int ego)
{
	std::vector<int> *setting = new std::vector<int>;
	std::set<int> neighbors;
	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		 iter.valid();
		 iter.next())
	{
		neighbors.insert(iter.actor());
	}
	for (IncidentTieIterator iter = pNetwork->inTies(ego);
		 iter.valid();
		 iter.next())
	{
		neighbors.insert(iter.actor());
	}
	neighbors.insert(ego);
	// when finished (not all done here) copy to the vector
	for (std::set<int>::const_iterator iter1 = neighbors.begin();
		 iter1 != neighbors.end(); iter1++)
	{
		setting->push_back(*iter1);
	}
	return setting;
}
}
