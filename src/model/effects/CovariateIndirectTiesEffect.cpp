/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateIndirectTiesEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateIndirectTiesEffect.
 *****************************************************************************/
#include <stdexcept>

#include "CovariateIndirectTiesEffect.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
CovariateIndirectTiesEffect::CovariateIndirectTiesEffect(
	const EffectInfo * pEffectInfo) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateIndirectTiesEffect::calculateContribution(int alter)
	const
{
	double change = 0;

	// If there are enough two-paths from the ego i to the alter j, then
	// we loose the distance 2 pair (i,j) by introducing the tie between
	// them.

	if (this->pTwoPathTable()->get(alter) != 0)
	{
		change -= this->value(alter);
	}

	// This variable is to simplify the later tests if a two-path through
	// the given alter makes a difference.

	int criticalTwoPathCount = 0;

	if (this->outTieExists(alter))
	{
		criticalTwoPathCount = 1;
	}

	// Consider each outgoing tie of the alter j.

	for (IncidentTieIterator iter = this->pNetwork()->outTies(alter);
		iter.valid();
		iter.next())
	{
		int h = iter.actor();

		// If h is not the ego i, there's no tie from i to h, and the
		// introduction or withdrawal of the tie (i,j) makes a difference
		// for the pair <i,h> to be a valid distance two pair,
		// then increment the contribution.

		if (h != this->ego() &&
			!this->outTieExists(h) &&
			this->pTwoPathTable()->get(h) == criticalTwoPathCount)
		{
			change += this->value(h);
		}
	}

	return change;
}


/**
 * This method is called at the start of the calculation of the statistic.
 */
void CovariateIndirectTiesEffect::initializeStatisticCalculation()
{
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();

	// A helper array of marks

	this->lmark = new int[n];

	for (int i = 0; i < n; i++)
	{
		this->lmark[i] = -1;
	}
}


/**
 * This method is called at the end of the calculation of the statistic.
 */
void CovariateIndirectTiesEffect::cleanupStatisticCalculation()
{
	delete[] this->lmark;
}


/**
 * Calculates the statistic corresponding to the given ego. The parameter
 * pNetwork is always the current network as there are no endowment effects
 * of this kind.
 */
double CovariateIndirectTiesEffect::egoStatistic(int ego,
	const Network * pNetwork)
{
	double statistic = 0;

	const Network * pStartMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period());
	const Network * pEndMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period() + 1);

	int i = ego;

	// Invariant: mark[h] = i if and only if a two-path from i
	// to h has been found.

	// Traverse all two-paths from i

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

			if (this->lmark[h] < i)
			{
				// The first two-path from i to h is found.

				this->lmark[h] = i;
				statistic += this->value(h);
			}
		}
	}

	// Okay, if there's a tie (i,h) then <i,h> cannot possibly be a
	// distance-two pair. Hence we iterate over outgoing ties (i,h) of i,
	// and if value(h) has been added to the statistic, we subtract it.

	for (IncidentTieIterator iter = pNetwork->outTies(i);
		iter.valid();
		iter.next())
	{
		int h = iter.actor();

		if (this->lmark[h] == i)
		{
			this->lmark[h] = -1;
			statistic -= this->value(h);
		}
	}

	// We do a similar fix for missing ties (i,h) at either end of
	// the period.

	for (IncidentTieIterator iter = pStartMissingNetwork->outTies(i);
		iter.valid();
		iter.next())
	{
		int h = iter.actor();

		if (this->lmark[h] == i)
		{
			this->lmark[h] = -1;
			statistic -= this->value(h);
		}
	}

	for (IncidentTieIterator iter = pEndMissingNetwork->outTies(i);
		iter.valid();
		iter.next())
	{
		int h = iter.actor();

		if (this->lmark[h] == i)
		{
			this->lmark[h] = -1;
			statistic -= this->value(h);
		}
	}

	// Ignore the trivial pair <i,i>.

	if (this->lmark[i] == i)
	{
		statistic -= this->value(i);
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double CovariateIndirectTiesEffect::endowmentStatistic(
	Network * pLostTieNetwork)
{
	throw logic_error(
		"CovariateIndirectTiesEffect: Endowment effect not supported");
}

}
