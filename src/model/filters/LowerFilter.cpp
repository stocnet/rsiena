/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: LowerFilter.cpp
 *
 * Description: This file contains the implementation of the class LowerFilter.
 *****************************************************************************/

#include "LowerFilter.h"
#include "network/IncidentTieIterator.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/NetworkChange.h"

namespace siena
{

/**
 * Constructs a new filter.
 */
LowerFilter::LowerFilter(const NetworkVariable * pOwnerVariable,
	const NetworkVariable * pOtherVariable) :
		NetworkDependentFilter(pOwnerVariable, pOtherVariable)
{
}


/**
 * Forbids tie changes between the given ego and some alters
 * by setting the permitted flag to false for these alters.
 */
void LowerFilter::filterPermittedChanges(int ego, bool * permitted)
{
	const Network * pNetwork1 = this->pVariable()->pNetwork();
	const Network * pNetwork2 = this->pOtherVariable()->pNetwork();

	// We shouldn't introduce a tie if it is not present in the other network.

	IncidentTieIterator iter1 = pNetwork1->outTies(ego);
	IncidentTieIterator iter2 = pNetwork2->outTies(ego);

	for (int j = 0; j < pNetwork1->m(); j++)
	{
		while (iter1.valid() && iter1.actor() < j)
		{
			iter1.next();
		}

		while (iter2.valid() && iter2.actor() < j)
		{
			iter2.next();
		}

		bool tieExists1 = iter1.valid() && iter1.actor() == j;
		bool tieExists2 = iter2.valid() && iter2.actor() == j;

		int noChange = ego;
		if (!this->pVariable()->oneModeNetwork())
		{
			noChange = pNetwork1->m();
		}
		if (!tieExists1 && !tieExists2 && j != noChange)
		{
			permitted[j] = false;
		}
	}
}


/**
 * Returns if applying the given ministep on the current state of the
 * network would be valid with respect to this filter.
 */
bool LowerFilter::validMiniStep(const NetworkChange * pMiniStep)
{
	const Network * pNetwork1 = this->pVariable()->pNetwork();
	const Network * pNetwork2 = this->pOtherVariable()->pNetwork();

	// We shouldn't introduce a tie if it is not present in the other network.

	int i = pMiniStep->ego();
	int j = pMiniStep->alter();

	return pNetwork1->tieValue(i, j) || pNetwork2->tieValue(i, j);
}

}
