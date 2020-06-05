/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AtLeastOneFilter.cpp
 *
 * Description: This file contains the implementation of the class
 * AtLeastOneFilter.
 *****************************************************************************/

#include "AtLeastOneFilter.h"
#include "network/IncidentTieIterator.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/NetworkChange.h"

namespace siena
{

/**
 * Constructs a new filter.
 */
AtLeastOneFilter::AtLeastOneFilter(const NetworkVariable * pOwnerVariable,
	const NetworkVariable * pOtherVariable) :
		NetworkDependentFilter(pOwnerVariable, pOtherVariable)
{
}


/**
 * Forbids tie changes between the given ego and some alters
 * by setting the permitted flag to false for these alters.
 */
void AtLeastOneFilter::filterPermittedChanges(int ego, bool * permitted)
{
	const Network * pNetwork1 = this->pVariable()->pNetwork();
	const Network * pNetwork2 = this->pOtherVariable()->pNetwork();

	// We shouldn't withdraw a tie if it is not present in the other network.

	IncidentTieIterator iter1 = pNetwork1->outTies(ego);
	IncidentTieIterator iter2 = pNetwork2->outTies(ego);

	while (iter1.valid())
	{
		while (iter2.valid() && iter2.actor() < iter1.actor())
		{
			iter2.next();
		}

		if (!iter2.valid() || iter2.actor() > iter1.actor())
		{
			permitted[iter1.actor()] = false;
		}

		iter1.next();
	}
}


/**
 * Returns if applying the given ministep on the current state of the
 * network would be valid with respect to this filter.
 */
bool AtLeastOneFilter::validMiniStep(const NetworkChange * pMiniStep)
{
	const Network * pNetwork1 = this->pVariable()->pNetwork();
	const Network * pNetwork2 = this->pOtherVariable()->pNetwork();

	// We shouldn't withdraw a tie if it is not present in the other network.

	int i = pMiniStep->ego();
	int j = pMiniStep->alter();

	return !pNetwork1->tieValue(i, j) || pNetwork2->tieValue(i, j);
}

}
