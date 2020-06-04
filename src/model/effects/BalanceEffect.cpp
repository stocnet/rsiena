/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BalanceEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BalanceEffect.
 *****************************************************************************/

#include <stdexcept>
#include "BalanceEffect.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "network/Network.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"
#include "model/tables/ConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
BalanceEffect::BalanceEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
	this->lbalanceMean = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void BalanceEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);

	const OneModeNetworkLongitudinalData * pNetworkData =
		dynamic_cast<const OneModeNetworkLongitudinalData *>(this->pData());

	if (pNetworkData)
	{
		this->lbalanceMean = pNetworkData->balanceMean();
	}
	else
	{
		throw logic_error("Data for one-mode network variable '" +
			this->pEffectInfo()->variableName() +
			"' expected.");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double BalanceEffect::calculateContribution(int alter) const
{
	// The formula from SIENA manual:
	// s_i(x) = sum_j x_{ij} sum_{h!=i,j} (b_0 - |x_{ih} - x_{jh}|)
	// Rearrange and obtain s_i(x) = A - B, where
	//    A = x_{i+} (n-2) b_0 and
	//    B = sum_j x_{ij} sum_{h!=i,j} |x_{ih} - x_{jh}|
	// A is easy to compute.
	// B equals the number of non-transitive two-paths starting at i
	// plus the number of out-stars <(i,j),(i,h)> with no tie from
	// j to h.

	// The contribution of a tie flip to A is (n-2) b_0.
	double a = (this->pNetwork()->n() - 2) * this->lbalanceMean;

	// These will be used later.

	int twoPathCount = this->pTwoPathTable()->get(alter);
	int inStarCount = this->pInStarTable()->get(alter);

	// First consider how the number of non-transitive two-paths would
	// change after introducing the tie from the ego i to the alter j.

	// x_{j+} - x_{ji} - IS_{ij} non-transitive two-paths would be created.

	double b =
		this->pNetwork()->outDegree(alter) - inStarCount;

	if (this->inTieExists(alter))
	{
		b--;
	}

	// However, TP_{ij} non-transitive two-paths would be lost as the
	// tie (i,j) would make them transitive.

	b -= twoPathCount;

	// The number of ties from the ego to actors other than the alter.

	int outDegree = this->pNetwork()->outDegree(this->ego());

	if (this->outTieExists(alter))
	{
		outDegree--;
	}

	// Now, consider the change in the number of out-stars <(i,j),(i,h)> with
	// no tie from j to h.

	// The introduction of the tie from the ego to the alter would create
	// x_{i+} - IS_{ij} new such out-stars, where the alter assumes the
	// role of j, and x_{i+} - TP_{ij} new out-stars, where the alter assumes
	// the role of h.

	b += 2 * outDegree - inStarCount - twoPathCount;

	return a - b;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double BalanceEffect::tieStatistic(int alter)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();
	const Network * pStartMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period());
	const Network * pEndMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period() + 1);

	// Initially all actors are valid

	this->lround++;
	int validActorCount = n;

	// Mark as invalid the actors h that have a missing tie from the ego or
	// the alter of the current tie at either end of the period.

	this->markInvalidActors(pStartMissingNetwork->outTies(this->ego()),
		validActorCount);
	this->markInvalidActors(pStartMissingNetwork->outTies(alter),
		validActorCount);
	this->markInvalidActors(pEndMissingNetwork->outTies(this->ego()),
		validActorCount);
	this->markInvalidActors(pEndMissingNetwork->outTies(alter),
		validActorCount);

	// Mark the ego and alter invalid as well.

	if (this->lflag[this->ego()] < this->lround)
	{
		this->lflag[this->ego()] = this->lround;
		validActorCount--;
	}

	if (this->lflag[alter] < this->lround)
	{
		this->lflag[alter] = this->lround;
		validActorCount--;
	}

	// Now we add the expression
	//   sum_h (b_0 - |x_{ih} - x_{jh}|)
	// to the statistic, where the sum extends over all valid actors h.

	// First add sum_h b_0
	statistic += validActorCount * this->lbalanceMean;

	// Now subtract sum_h |x_{ih} - x_{jh}|

	IncidentTieIterator egoIter = pNetwork->outTies(this->ego());
	IncidentTieIterator alterIter = pNetwork->outTies(alter);

	while (egoIter.valid() || alterIter.valid())
	{
		if (egoIter.valid() &&
			(!alterIter.valid() || egoIter.actor() < alterIter.actor()))
		{
			if (this->lflag[egoIter.actor()] < this->lround)
			{
				statistic--;
			}

			egoIter.next();
		}
		else if (alterIter.valid() &&
			(!egoIter.valid() || alterIter.actor() < egoIter.actor()))
		{
			if (this->lflag[alterIter.actor()] < this->lround)
			{
				statistic--;
			}

			alterIter.next();
		}
		else
		{
			egoIter.next();
			alterIter.next();
		}
	}

	return statistic;
}


/**
 * This method is called at the start of the calculation of the statistic.
 */
void BalanceEffect::initializeStatisticCalculation()
{
	int n = this->pNetwork()->n();

	this->lflag = new int[n];
	this->lround = 0;

	for (int i = 0; i < n; i++)
	{
		this->lflag[i] = 0;
	}
}


/**
 * This method is called at the end of the calculation of the statistic.
 */
void BalanceEffect::cleanupStatisticCalculation()
{
	delete[] this->lflag;
}


/**
 * This method marks as invalid all actors that are iterated over by the given
 * iterator. The fact that an actor i is invalid is represented by setting
 * lflag[i] = lround. It is assumed that lflag[i] <= lround for all actors,
 * meaning that lflag[i] < lround holds for valid actors. The variable
 * validActorCount keeps track of the still valid actors.
 */
void BalanceEffect::markInvalidActors(IncidentTieIterator iter,
	int & validActorCount)
{
	while (iter.valid())
	{
		if (this->lflag[iter.actor()] < this->lround)
		{
			this->lflag[iter.actor()] = this->lround;
			validActorCount--;
		}

		iter.next();
	}
}

}
