/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BetweennessEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BetweennessEffect.
 *****************************************************************************/

#include "BetweennessEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
BetweennessEffect::BetweennessEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double BetweennessEffect::calculateContribution(int alter) const
{
	// Get the number of actors h other than the alter j that have a tie
	// to the ego i.

	int inDegree = this->pNetwork()->inDegree(this->ego());

	if (this->inTieExists(alter))
	{
		inDegree--;
	}

	// Now, for each of these actors h, the introduction of the tie (i,j)
	// creates a new non-transitive two-path through i unless there's a tie
	// (h,j), in which case <(h,i),(h,j)> is an out-star.

	return inDegree - this->pOutStarTable()->get(alter);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double BetweennessEffect::tieStatistic(int alter)
{
	int statistic = 0;
	const Network * pNetwork = this->pNetwork();
	const Network * pStartMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period());
	const Network * pEndMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period() + 1);

	// Calculate the number of non-transitive two-paths through
	// the ego i terminating at the alter j.

	int i = this->ego();
	int j = alter;

	// For the start, assume that each in-neighbor of i gives one such
	// a two-path.

	statistic += pNetwork->inDegree(i);

	// The actor j itself doesn't count

	if (this->lmark[j] >= this->lbaseMark)
	{
		statistic--;
	}

	// There are three cases where an actor h shoudln't contribute:
	// - There is a tie (h,j).
	// - The tie (h,j) is missing at the start of the period.
	// - The tie (h,j) is missing at the end of the period.
	// We should decrement the statistic for each of these actors h,
	// but we shouldn't do that more than once per actor, so we use
	// another mark.

	this->lcurrentMark++;

	// Ties (h,j)

	for (IncidentTieIterator iterH = pNetwork->inTies(j);
		iterH.valid();
		iterH.next())
	{
		int h = iterH.actor();

		if (this->lmark[h] >= this->lbaseMark &&
			this->lmark[h] < this->lcurrentMark)
		{
			statistic--;
			this->lmark[h] = this->lcurrentMark;
		}
	}

	// Missing ties (h,j) at the start of the period

	for (IncidentTieIterator iterH = pStartMissingNetwork->inTies(j);
		iterH.valid();
		iterH.next())
	{
		int h = iterH.actor();

		if (this->lmark[h] >= this->lbaseMark &&
			this->lmark[h] < this->lcurrentMark)
		{
			statistic--;
			this->lmark[h] = this->lcurrentMark;
		}
	}

	// Missing ties (h,j) at the end of the period

	for (IncidentTieIterator iterH = pEndMissingNetwork->inTies(j);
		iterH.valid();
		iterH.next())
	{
		int h = iterH.actor();

		if (this->lmark[h] >= this->lbaseMark &&
			this->lmark[h] < this->lcurrentMark)
		{
			statistic--;
			this->lmark[h] = this->lcurrentMark;
		}
	}

	return statistic;
}


/**
 * This method is called at the start of the calculation of the statistic.
 */
void BetweennessEffect::initializeStatisticCalculation()
{
	int n = this->pNetwork()->n();

	// A helper array of marks

	this->lmark = new int[n];
	this->lcurrentMark = 0;

	for (int i = 0; i < n; i++)
	{
		this->lmark[i] = 0;
	}
}


/**
 * This method is called at the end of the calculation of the statistic.
 */
void BetweennessEffect::cleanupStatisticCalculation()
{
	delete[] this->lmark;
}


/**
 * This method is called right before summing up the contributions of the
 * outgoing ties of the given ego in the calculation of the statistic.
 */
void BetweennessEffect::onNextEgo(int ego)
{
	const Network * pNetwork = this->pNetwork();

	// Mark the in-neighbors of the ego i by ensuring the following assertion:
	// mark[h] >= currentMark <==> the tie (h,i) exists

	this->lcurrentMark++;

	for (IncidentTieIterator iter = pNetwork->inTies(ego);
		iter.valid();
		iter.next())
	{
		this->lmark[iter.actor()] = this->lcurrentMark;
	}

	// Remember the current mark such that mark[h] >= baseMark if and only
	// if there's a tie from h to i.

	this->lbaseMark = this->lcurrentMark;
}

}
