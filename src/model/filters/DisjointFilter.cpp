/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DisjointFilter.cpp
 *
 * Description: This file contains the implementation of the class
 * DisjointFilter.
 *****************************************************************************/

#include "DisjointFilter.h"
#include "network/IncidentTieIterator.h"
#include "network/iterators/UnionTieIterator.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/NetworkChange.h"

namespace siena
{

/**
 * Constructs a new filter.
 */
DisjointFilter::DisjointFilter(const NetworkVariable * pOwnerVariable,
	const NetworkVariable * pOtherVariable) :
		NetworkDependentFilter(pOwnerVariable, pOtherVariable)
{
	const NetworkVariable * pNetworkVariable1 = this->pVariable();
	const NetworkVariable * pNetworkVariable2 = this->pOtherVariable();
	this->lsymm = ((pNetworkVariable1->symmetric()) && (!pNetworkVariable2->symmetric()));
//  The extra steps for a symmetric owner network are superfluous
//  if the other network is also symmetric.
}


/**
 * Forbids tie changes between the given ego and some alters
 * by setting the permitted flag to false for these alters.
 */
void DisjointFilter::filterPermittedChanges(int ego, bool * permitted)
{
	const Network * pNetwork1 = this->pVariable()->pNetwork();
	const Network * pNetwork2 = this->pOtherVariable()->pNetwork();

	// We shouldn't introduce a tie if it is present in the other network.

	IncidentTieIterator iter1 = pNetwork1->outTies(ego);
	IncidentTieIterator iter2 = pNetwork2->outTies(ego);

	while (iter2.valid())
	{
		while (iter1.valid() && iter1.actor() < iter2.actor())
		{
			iter1.next();
		}

		if (!iter1.valid() || iter1.actor() > iter2.actor())
		{
			permitted[iter2.actor()] = false;
		}

		iter2.next();
	}

	if (this->lsymm)
	{
		IncidentTieIterator iter1 = pNetwork1->outTies(ego);
		IncidentTieIterator iter2 = pNetwork2->inTies(ego, "djf");
		while (iter2.valid())
		{
			while (iter1.valid() && iter1.actor() < iter2.actor())
			{
				iter1.next();
			}

			if (!iter1.valid() || iter1.actor() > iter2.actor())
			{
				permitted[iter2.actor()] = false;
			}

			iter2.next();
		}
	}
}

/**
 * Returns if applying the given ministep on the current state of the
 * network would be valid with respect to this filter.
 */
bool DisjointFilter::validMiniStep(const NetworkChange * pMiniStep)
{
	const Network * pNetwork1 = this->pVariable()->pNetwork();
	const Network * pNetwork2 = this->pOtherVariable()->pNetwork();

	// We shouldn't introduce a tie if it is present in the other network.

	int i = pMiniStep->ego();
	int j = pMiniStep->alter();

	if  (this->lsymm)
	{
		return (pNetwork1->tieValue(i, j) ||
				!(pNetwork2->tieValue(i, j) || pNetwork2->tieValue(j, i)));
	}
	else
	{
		return pNetwork1->tieValue(i, j) || !pNetwork2->tieValue(i, j);
	}
}

}
