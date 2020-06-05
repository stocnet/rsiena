/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedTwoPathTable.cpp
 *
 * Description: This file contains the implementation of the MixedTwoPathTable
 * class.
 *****************************************************************************/
#include <R_ext/Print.h>

#include "MixedTwoPathTable.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "network/OneModeNetwork.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Creates a new table for storing two-paths with the specified directions for
 * the first and the second tie.
 */
MixedTwoPathTable::MixedTwoPathTable(TwoNetworkCache * pOwner,
	Direction firstStepDirection,
	Direction secondStepDirection) : MixedEgocentricConfigurationTable(pOwner)
{
	this->lfirstStepDirection = firstStepDirection;
	this->lsecondStepDirection = secondStepDirection;
}


// ----------------------------------------------------------------------------
// Section: ConfigurationTable implementation
// ----------------------------------------------------------------------------

/**
 * Calculates the number of generalized two-paths between the ego and all
 * other actors.
 */
void MixedTwoPathTable::calculate()
{
	// Reset the counters to zeroes
	this->reset();
	if (this->lfirstStepDirection == FORWARD)
	{
		this->performFirstStep(
				this->pFirstNetwork()->outTies(this->ego()));
	}
	else if (this->lfirstStepDirection == BACKWARD)
	{
		this->performFirstStep(
				this->pFirstNetwork()->inTies(this->ego()));
	}
	else
	{
		const OneModeNetwork * pOneModeFirstNetwork =
			dynamic_cast<const OneModeNetwork *>(this->pFirstNetwork());
		this->performFirstStep(
				pOneModeFirstNetwork->reciprocatedTies(this->ego()));
	}
}


/**
 * Performs the first step by iterating over the actors of the given
 * iterator and invoking the method for the second step.
 */
template<class Iterator>
void MixedTwoPathTable::performFirstStep(Iterator iter)
{
	// TODO: Using templates here is a bad design. It's because
	// IncidentTieIterator and CommonNeighborIterator have no
	// common base class.

	// Try out all possible first steps

	while (iter.valid())
	{
		int middleActor = iter.actor();
		iter.next();
		// Choose the right iterator for the second step
		if (this->lsecondStepDirection == FORWARD)
		{
			this->performSecondStep(
					this->pSecondNetwork()->outTies(middleActor));
		}
		else if (this->lsecondStepDirection == BACKWARD)
		{
			this->performSecondStep(
					this->pSecondNetwork()->inTies(middleActor));
		}
		else
		{
			const OneModeNetwork * pOneModeSecondNetwork =
				dynamic_cast<const OneModeNetwork *>(this->pSecondNetwork());
			this->performSecondStep(
					pOneModeSecondNetwork->reciprocatedTies(middleActor));
		}
	}
}


/**
 * Performs the second step by iterating over the actors of the given
 * iterator and incrementing their values that are stored in this table.
 */
template<class Iterator>
void MixedTwoPathTable::performSecondStep(Iterator iter)
{
	while (iter.valid())
	{
		this->ltable[iter.actor()]++;
		iter.next();
	}
}

}
