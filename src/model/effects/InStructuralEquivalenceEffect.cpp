/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InStructuralEquivalenceEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * InStructuralEquivalenceEffect.
 *****************************************************************************/

#include <stdexcept>
#include "InStructuralEquivalenceEffect.h"
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
InStructuralEquivalenceEffect::InStructuralEquivalenceEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
/*
  This value needs to be changed later if we wish to center the effects;
  later on commented out by //
	this->linStructEqMean = 0;
*/
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void InStructuralEquivalenceEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);

	const OneModeNetworkLongitudinalData * pNetworkData =
		dynamic_cast<const OneModeNetworkLongitudinalData *>(this->pData());


	if (pNetworkData)
	{
		this->linStructEqMean = pNetworkData->structuralMean();
	}
	else
	{
		throw std::logic_error("Data for one-mode network variable '" +
			this->pEffectInfo()->variableName() + "' expected.");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double InStructuralEquivalenceEffect::calculateContribution(int alter) const
{   return
      (this->pNetwork()->n() - 2)*linStructEqMean
       - this->pNetwork()->inDegree(this->ego() )
       - this->pNetwork()->inDegree(alter)
       + this->inTieExists(alter)
       + this->outTieExists(alter)
       + 2 * this->pOutStarTable()->get(alter);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double InStructuralEquivalenceEffect::tieStatistic(int alter)
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

	this->markInvalidActors(pStartMissingNetwork->inTies(this->ego()),
		validActorCount);
	this->markInvalidActors(pStartMissingNetwork->inTies(alter),
		validActorCount);
	this->markInvalidActors(pEndMissingNetwork->inTies(this->ego()),
		validActorCount);
	this->markInvalidActors(pEndMissingNetwork->inTies(alter),
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
	//   sum_h (1 - |x_{hi} - x_{hj}|)
	// to the statistic, wthe sum extends over all valid actors h.

	// First add sum_h b_0
   	statistic += validActorCount * this->linStructEqMean;

	// Now subtract sum_h |x_{hi} - x_{hj}|

	IncidentTieIterator egoIter = pNetwork->inTies(this->ego());
	IncidentTieIterator alterIter = pNetwork->inTies(alter);

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
void InStructuralEquivalenceEffect::initializeStatisticCalculation()
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
void InStructuralEquivalenceEffect::cleanupStatisticCalculation()
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
void InStructuralEquivalenceEffect::markInvalidActors(IncidentTieIterator iter,
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
