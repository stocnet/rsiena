/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Chain.cpp
 *
 * Description: This file contains the implementation of the class Chain.
 *****************************************************************************/

#include <vector>
#include <stdexcept>
#include <string>

#include <Rinternals.h>

#include "Chain.h"
#include "utils/Utils.h"
#include "utils/Random.h"
#include "model/State.h"
#include "model/ml/MiniStep.h"
#include "model/ml/BehaviorChange.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/MLSimulation.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/Data.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/NetworkLongitudinalData.h"
#include "model/variables/DependentVariable.h"
#include <R_ext/Error.h>
#include <Rinternals.h>

using namespace std;

namespace siena
{
SEXP getChainDF(const Chain& chain, bool sort);
SEXP getMiniStepDF(const MiniStep& miniStep);

// ----------------------------------------------------------------------------
// Section: Constructors, destructors, initializers
// ----------------------------------------------------------------------------

/**
 * Creates an empty chain.
 */
Chain::Chain(Data * pData)
{
	this->lpFirst = new MiniStep(0, 0);
	this->lpLast = new MiniStep(0, 0);

	this->lpFirst->pNext(this->lpLast);
	this->lpLast->pPrevious(this->lpFirst);

	this->resetOrderingKeys();

	this->lpData = pData;
	this->lperiod = -1;

	this->lpInitialState = 0;

	this->lminiSteps.push_back(this->lpLast);
	this->lpLast->index(0);

	this->lmu = 0;
	this->lsigma2 = 0;
	this->lfinalReciprocalRate = 0;
}


/**
 * Deallocates this chain.
 */
Chain::~Chain()
{
	this->clear();

	delete this->lpFirst;
	delete this->lpLast;
	this->lpFirst = 0;
	this->lpLast = 0;
	this->lpData = 0;

	this->lminiSteps.clear();
	delete this->lpInitialState;
	this->lpInitialState = 0;
	deallocateVector(this->linitialStateDifferences);
	deallocateVector(this->lendStateDifferences);
}


/**
 * Stores the initial state as of the beginning of the current period.
 * The parameter <code>copyValues</code> indicates if the values of the
 * Data object should be copied or simply referenced.
 */
void Chain::setupInitialState(bool copyValues)
{
	delete this->lpInitialState;
	this->lpInitialState =
		new State(this->lpData, this->lperiod, copyValues);
}


// ----------------------------------------------------------------------------
// Section: Chain modifications
// ----------------------------------------------------------------------------

/**
 * Removes all ministeps from this chain except for the dummy ministeps
 * at the ends of the chain. NB Does not alter the initialState
 */
void Chain::clear()
{
	MiniStep * pMiniStep = this->lpFirst->pNext();

	while (pMiniStep != this->lpLast)
	{
		MiniStep * pNextMiniStep = pMiniStep->pNext();
		delete pMiniStep;
		pMiniStep = pNextMiniStep;
	}

	this->lpFirst->pNext(this->lpLast);
	this->lpLast->pPrevious(this->lpFirst);

	this->lminiSteps.clear();
	this->lminiSteps.push_back(this->lpLast);
	this->lpLast->index(0);

	this->ldiagonalMiniSteps.clear();
	this->lccpMiniSteps.clear();
	this->lmissingNetworkMiniSteps.clear();
	this->lmissingBehaviorMiniSteps.clear();
	this->lfirstMiniStepPerOption.clear();

	this->lmu = 0;
	this->lsigma2 = 0;
	this->lfinalReciprocalRate = 0;
}


/**
 * Inserts a new ministep before a ministep that already belongs to this chain.
 */
void Chain::insertBefore(MiniStep * pNewMiniStep, MiniStep * pExistingMiniStep)
{
	MiniStep * pPreviousMiniStep = pExistingMiniStep->pPrevious();

	pNewMiniStep->pChain(this);

	// Update the pointers to next and previous ministeps.

	pNewMiniStep->pPrevious(pPreviousMiniStep);
	pPreviousMiniStep->pNext(pNewMiniStep);

	pNewMiniStep->pNext(pExistingMiniStep);
	pExistingMiniStep->pPrevious(pNewMiniStep);

	// Set the ordering key.

	double leftKey = pPreviousMiniStep->orderingKey();
	double rightKey = pExistingMiniStep->orderingKey();
	double newKey = (leftKey + rightKey) / 2;

	if ((-leftKey + newKey) > 1e-10 && (- newKey +rightKey) > 1e-10)
	{
		pNewMiniStep->orderingKey(newKey);
	}
	else
	{
		// Bad. We have exhausted the double precision.
		// However, using a linear time effort, we can restore a correct order.

		this->resetOrderingKeys();
	}

	// Update the array of diagonal ministeps

	if (pNewMiniStep->diagonal())
	{
		this->ldiagonalMiniSteps.push_back(pNewMiniStep);
		pNewMiniStep->diagonalIndex(this->ldiagonalMiniSteps.size() - 1);
	}

	if (pNewMiniStep->missingEnd(this->lperiod))
//	if (pNewMiniStep->missing(this->lperiod))
	{
		if (pNewMiniStep->networkMiniStep())
		{
			pNewMiniStep->missingIndex(this->lmissingNetworkMiniSteps.size());
			this->lmissingNetworkMiniSteps.push_back(pNewMiniStep);
		}
		else
		{
			pNewMiniStep->missingIndex(this->lmissingBehaviorMiniSteps.size());
				this->lmissingBehaviorMiniSteps.push_back(pNewMiniStep);
		}
	}

	// Update the array of all ministeps

	this->lminiSteps.push_back(pNewMiniStep);
	pNewMiniStep->index(this->lminiSteps.size() - 1);

	this->updateSameOptionPointersOnInsert(pNewMiniStep);

	// Update the vector of CCPs for each ministep whose CCP status
	// might have changed.

	this->updateCCPs(pNewMiniStep->pPrevious());
	this->updateCCPs(pNewMiniStep->pPreviousWithSameOption());
	this->updateCCPs(pNewMiniStep);

	// Update some variables that are constantly maintained

	double rr = pNewMiniStep->reciprocalRate();
	this->lmu += rr;
	this->lsigma2 += rr * rr;
	// Rprintf("insert step %f %f %f\n", rr, this->lmu, this->lsigma2);
}


/**
 * Removes the given ministep from this chain without deleting the ministep.
 */
void Chain::remove(MiniStep * pMiniStep)
{
	MiniStep * pPrevious = pMiniStep->pPrevious();
	MiniStep * pPreviousWithSameOption = pMiniStep->pPreviousWithSameOption();

	// Updates pointers to the next and previous ministep.

	pPrevious->pNext(pMiniStep->pNext());
	pMiniStep->pNext()->pPrevious(pPrevious);
	pMiniStep->pNext(0);
	pMiniStep->pPrevious(0);

	// Next and previous pointers to ministeps of the same option

	if (this->lfirstMiniStepPerOption[*pMiniStep->pOption()] ==
		pMiniStep)
	{
		this->lfirstMiniStepPerOption[*pMiniStep->pOption()] =
			pMiniStep->pNextWithSameOption();
	}

	if (pMiniStep->pNextWithSameOption())
	{
		pMiniStep->pNextWithSameOption()->pPreviousWithSameOption(
			pPreviousWithSameOption);
	}

	if (pPreviousWithSameOption)
	{
		pPreviousWithSameOption->pNextWithSameOption(
			pMiniStep->pNextWithSameOption());
	}

	pMiniStep->pPreviousWithSameOption(0);
	pMiniStep->pNextWithSameOption(0);

	// Update the vector of CCPs

	this->updateCCPs(pPrevious);
	this->updateCCPs(pPreviousWithSameOption);
	this->updateCCPs(pMiniStep);

	// Update the vector of diagonal ministeps

	if (pMiniStep->diagonal())
	{
		MiniStep * pLastMiniStep =
			this->ldiagonalMiniSteps[this->ldiagonalMiniSteps.size() - 1];
		this->ldiagonalMiniSteps[pMiniStep->diagonalIndex()] = pLastMiniStep;
		pLastMiniStep->diagonalIndex(pMiniStep->diagonalIndex());
		this->ldiagonalMiniSteps.pop_back();
		pMiniStep->diagonalIndex(-1);
	}

	// Update the vectors of ministeps for missing options

	if (pMiniStep->missingEnd(this->lperiod))
//	if (pMiniStep->missing(this->lperiod))
	{
		vector<MiniStep *> * pVector = &this->lmissingNetworkMiniSteps;

		if (pMiniStep->behaviorMiniStep())
		{
			pVector = &this->lmissingBehaviorMiniSteps;
		}

		MiniStep * pLast = (*pVector)[pVector->size() - 1];
		(*pVector)[pMiniStep->missingIndex()] = pLast;
		pLast->missingIndex(pMiniStep->missingIndex());
		pVector->pop_back();
		pMiniStep->missingIndex(-1);
	}

	MiniStep * pLastMiniStep = this->lminiSteps[this->lminiSteps.size() - 1];
	this->lminiSteps[pMiniStep->index()] = pLastMiniStep;
	pLastMiniStep->index(pMiniStep->index());
	this->lminiSteps.pop_back();
	pMiniStep->index(-1);

	double rr = pMiniStep->reciprocalRate();
	this->lmu -= rr;
	this->lsigma2 -= rr * rr;

	pMiniStep->pChain(0);
}


/**
 * Generates a random chain connecting the start and end observations of the
 * given data object for the given period. The chain is simple in the sense
 * that no two ministeps cancel each other out.
 */
void Chain::connect(int period, MLSimulation * pMLSimulation)
{
	this->clear();
	this->lperiod = period;
	vector<MiniStep *> miniSteps;

	// Create the required ministeps

	for (unsigned variableIndex = 0;
		variableIndex < this->lpData->rDependentVariableData().size();
		variableIndex++)
	{
		LongitudinalData * pVariableData =
			this->lpData->rDependentVariableData()[variableIndex];
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(pVariableData);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(pVariableData);

		if (pNetworkData)
		{
			const Network * pNetwork1 = pNetworkData->pNetwork(period);
			const Network * pNetwork2 = pNetworkData->pNetwork(period + 1);

			for (int i = 0; i < pNetwork1->n(); i++)
			{
				IncidentTieIterator iter1 = pNetwork1->outTies(i);
				IncidentTieIterator iter2 = pNetwork2->outTies(i);

				while (iter1.valid() || iter2.valid())
				{
					if (iter1.valid() &&
						(!iter2.valid() || iter1.actor() < iter2.actor()))
					{
						if (!pNetworkData->structural(i, iter1.actor(),
								period)
							//	|| !pNetworkData->structural(i, iter1.actor(),
							//	period + 1)
							)
						{
							miniSteps.push_back(
								new NetworkChange(pNetworkData,
									i,
									iter1.actor(), false));
							iter1.next();
						}
						else
						{
							// create step in structural subchain?
						}
					}
					else if (iter2.valid() &&
						(!iter1.valid() || iter2.actor() < iter1.actor()))
					{
						if (!pNetworkData->structural(i, iter2.actor(),
								period)
							//	|| !pNetworkData->structural(i, iter2.actor(),
							// period + 1)
							)
						{
							miniSteps.push_back(
								new NetworkChange(pNetworkData,
									i,
									iter2.actor(), false));
							iter2.next();
						}
						else
						{
							// create step in structural subchain?
						}
					}
					else
					{
						iter1.next();
						iter2.next();
					}
				}
			}
		}
		else if (pBehaviorData)
		{
			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				int delta = pBehaviorData->value(period + 1, i) -
					pBehaviorData->value(period, i);
				int singleChange = 1;

				if (delta < 0)
				{
					delta = -delta;
					singleChange = -1;
				}

				for (int j = 0; j < delta; j++)
				{
					if (!pBehaviorData->structural(period, j)
						//|| !pBehaviorData->structural(period, j + 1)
						)

					{
						miniSteps.push_back(
							new BehaviorChange(pBehaviorData,
								i,
								singleChange));
					}
					else
					{
						// create step in structural subchain?
					}
				}
			}
		}
	}

	// If we have no constraints we can go ahead and add to the chain in a
	// random order. If we do have contraints we need to make sure our chain
	// is valid. For now we have two distinct loops:

	if (this->lpData->rNetworkConstraints().size()  == 0)
	{
		// Randomize the ministeps

		for (unsigned i = 1; i < miniSteps.size(); i++)
		{
			int j = nextInt(i + 1);
			MiniStep * pTempMiniStep = miniSteps[i];
			miniSteps[i] = miniSteps[j];
			miniSteps[j] = pTempMiniStep;
		}

		// And finally add the ministeps to this chain

		for (unsigned i = 0; i < miniSteps.size(); i++)
		{
			this->insertBefore(miniSteps[i], this->lpLast);
		}
	}
	else
	{
		unsigned count = 0;
		vector<MiniStep *> remainingMiniSteps;
		while(miniSteps.size() > 0 &&
			count < this->lpData->rDependentVariableData().size())
		{
			count++;
			remainingMiniSteps.clear();
			for (unsigned i = 0; i < miniSteps.size(); i++)
			{
				// first try random insert
				//const NetworkChange * pNetworkChange =
				//	dynamic_cast<const NetworkChange *>(miniSteps[i]);
				MiniStep * pMiniStep =
					this->randomMiniStep(this->lpFirst->pNext(),
						this->lpLast);
				bool valid = false;
				if (miniSteps[i]->behaviorMiniStep())
				{
					valid = true;
				}
				else
				{
					// get current state at this place
					pMLSimulation->initialize(this->lperiod);
					pMLSimulation->executeMiniSteps(this->lpFirst->pNext(),
						pMiniStep);
					// see if valid here
					DependentVariable * pVariable =
						pMLSimulation->rVariables()[miniSteps[i]->variableId()];
					//	Rf_PrintValue(getMiniStepDF(*miniSteps[i]));

					if (!pVariable->validMiniStep(miniSteps[i]))
					{
						//		Rprintf("first inval\n");
						// no go, so try to find somewhere else to put it
						MiniStep * pLastMiniStep =
							this->pLastMiniStepForLink(miniSteps[i]);
						if (pLastMiniStep != this->lpFirst)
						{
							//		Rprintf("lastl\n");
							//	Rf_PrintValue(getMiniStepDF(*pLastMiniStep));
							pMiniStep =
								this->randomMiniStep(pLastMiniStep->pNext(),
									this->lpLast);
							//	Rf_PrintValue(getMiniStepDF(*pMiniStep));
							pMLSimulation->initialize(this->lperiod);
							pMLSimulation->executeMiniSteps(this->lpFirst->
								pNext(), pMiniStep);
							// see if valid here
							if (!pVariable->validMiniStep(miniSteps[i]))
							{
								//		Rprintf("second inval\n");
								// no go, so try to find somewhere else to put
								// it
								MiniStep * pFirstMiniStep =
									this->pFirstMiniStepForLink(miniSteps[i]);
								//	if (pFirstMiniStep != this->lpFirst)
								// {
								//	Rf_PrintValue(getMiniStepDF(*pFirstMiniStep));
								pMiniStep =
									this->randomMiniStep(this->
										lpFirst->pNext(),
										pFirstMiniStep);
								//	Rf_PrintValue(getMiniStepDF(*pMiniStep));
								pMLSimulation->initialize(this->lperiod);
								pMLSimulation->
									executeMiniSteps(pFirstMiniStep,
										pMiniStep);
								// see if valid here
								if (pVariable->validMiniStep(miniSteps[i]))
								{
								//	Rprintf("third true val\n");
									valid = true;
								}
								//	}
							}
							else
							{
								valid = true;
							}
						}
					}
					else
					{
						valid = true;
					}
				}
				if (valid)
				{
					this->insertBefore(miniSteps[i], pMiniStep);
				}
				else
				{
					remainingMiniSteps.push_back(miniSteps[i]);
				}
			}
			miniSteps = remainingMiniSteps;
		}
		//	Rf_PrintValue(getChainDF(*this, false));
		//	Rprintf("****** count %d\n", count);
		if (miniSteps.size() > 0)
		{
			for (unsigned i = 0; i < miniSteps.size(); i++)
			{
				Rf_PrintValue(getMiniStepDF(*miniSteps[i]));
			}

			Rf_error("Cannot create minimal chain due to constraints");
		}
	}
}


// ----------------------------------------------------------------------------
// Section: Various updates
// ----------------------------------------------------------------------------

/**
 * Sets the value of the period of the chain.
 */
void Chain::period(int value)
{
	this->lperiod = value;
}

/**
 * Performs the necessary updates when the reciprocal rate of the
 * given ministep changes to the given value.
 */
void Chain::onReciprocalRateChange(const MiniStep * pMiniStep, double newValue)
{
	double oldValue = pMiniStep->reciprocalRate();

	this->lmu -= oldValue;
	this->lsigma2 -= oldValue * oldValue;

	this->lmu += newValue;
	this->lsigma2 += newValue * newValue;
}


/**
 * Changes the initial state according to the given ministep.
 */
void Chain::changeInitialState(const MiniStep * pMiniStep)
{
	//Rprintf("%d change\n",this->lperiod);
//	Rf_PrintValue(getMiniStepDF(*pMiniStep));
	if (pMiniStep->networkMiniStep())
	{
		const NetworkChange * pNetworkChange =
			dynamic_cast<const NetworkChange *>(pMiniStep);

		// Okay, this is a bad trick and indicates an imperfect design.
		// We cast away the constness because we need to change the network.
		// But this method is called only in situation when the initial
		// state is actually a copy of the observed values and not simply
		// references to networks in the Data object, so we don't destroy
		// the observed Data object.

		Network * pNetwork =
			(Network *) this->lpInitialState->pNetwork(
				pNetworkChange->variableName());
		int ego = pNetworkChange->ego();
		int alter = pNetworkChange->alter();
		pNetwork->setTieValue(ego, alter, 1 - pNetwork->tieValue(ego, alter));
	}
	else
	{
		const BehaviorChange * pBehaviorChange =
			dynamic_cast<const BehaviorChange *>(pMiniStep);

		// Same misdesign here

		int * values =
			(int *) this->lpInitialState->behaviorValues(
				pBehaviorChange->variableName());
		values[pBehaviorChange->ego()] += pBehaviorChange->difference();
	}
}

/**
 *  Adds a single ministep to the vector of initialStateDifferences
 */
void Chain::addInitialStateDifference(MiniStep * pMiniStep)
{
	this->linitialStateDifferences.push_back(pMiniStep);
}

/**
 *  Updates the initialState of the chain with the ministeps version of the
 *  differences.
 */

void Chain::recreateInitialState()
{
//	Rprintf("ff %d %d\n", this->linitialStateDifferences.size(), this->lperiod);
	for (unsigned i= 0; i < this->linitialStateDifferences.size(); i++)
	{
//	Rprintf("ff\n");
	this->changeInitialState(this->linitialStateDifferences[i]);
	}
}
/**
 * Generates a chain connecting the original and simulated initial
 * observations of the given data object for the given period.
 * The chain is simple in the sense that no two ministeps cancel each other out.
 */
void Chain::createInitialStateDifferences()
{
	deallocateVector(this->linitialStateDifferences);
	//Rprintf("******** %d create initial\n",this->linitialStateDifferences.size() );
	Data * pData = this->lpData;
	State * initialState = this->lpInitialState;
	int period = this->lperiod;

	// Create the required ministeps

	for (unsigned variableIndex = 0;
		variableIndex < pData->rDependentVariableData().size();
		variableIndex++)
	{
		LongitudinalData * pVariableData =
			pData->rDependentVariableData()[variableIndex];
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(pVariableData);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(pVariableData);

		if (pNetworkData)
		{
			const Network * pNetwork1 = pNetworkData->pNetwork(period);
			const Network * pNetwork2 =
				initialState->pNetwork(pNetworkData->name());

			for (int i = 0; i < pNetwork1->n(); i++)
			{
				IncidentTieIterator iter1 = pNetwork1->outTies(i);
				IncidentTieIterator iter2 = pNetwork2->outTies(i);

				while (iter1.valid() || iter2.valid())
				{
					if (iter1.valid() &&
						(!iter2.valid() || iter1.actor() < iter2.actor()))
					{
						if (!pNetworkData->structural(i, iter1.actor(),
								period)
							//	|| !pNetworkData->structural(i, iter1.actor(),
							//	period + 1)
							)
						{
							NetworkChange * pMiniStep =
								new NetworkChange(pNetworkData,
									i,
									iter1.actor(), false);
							//		Rf_PrintValue(getMiniStepDF(*pMiniStep));
							this->linitialStateDifferences.
								push_back(pMiniStep);

							iter1.next();
							//	Rf_PrintValue(getMiniStepDF(
							//		*this->linitialStateDifferences.back()));
						}
						else
						{
							// create step in structural subchain?
						}
					}
					else if (iter2.valid() &&
						(!iter1.valid() || iter2.actor() < iter1.actor()))
					{
						if (!pNetworkData->structural(i, iter2.actor(),
								period)
							//	|| !pNetworkData->structural(i, iter2.actor(),
							// period + 1)
							)
						{
							this->linitialStateDifferences.push_back(
								new NetworkChange(pNetworkData,
									i,
									iter2.actor(), false));
							iter2.next();
							//Rf_PrintValue(getMiniStepDF(
							//		*this->linitialStateDifferences.back()));
						}
						else
						{
							// create step in structural subchain?
						}
					}
					else
					{
						iter1.next();
						iter2.next();
					}
				}
			}
		}
		else if (pBehaviorData)
		{
			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				int delta = initialState->
					behaviorValues(pBehaviorData->name())[i]
					- pBehaviorData->value(period, i);
				int singleChange = 1;

				if (delta < 0)
				{
					delta = -delta;
					singleChange = -1;
				}

				for (int j = 0; j < delta; j++)
				{
					if (!pBehaviorData->structural(period, j)
						//|| !pBehaviorData->structural(period, j + 1)
						)

					{
						this->linitialStateDifferences.push_back(
							new BehaviorChange(pBehaviorData,
								i,
								singleChange));
						//	Rprintf(" %d %d in beh\n", i, singleChange);
					}
					else
					{
						// create step in structural subchain?
					}
				}
			}
		}
	}
//	Rprintf("xx ********%d %d end create initial diff\n", this->linitialStateDifferences.size(), period);
}
/**
 *  Adds a single ministep to the vector of endStateDifferences
 */
void Chain::addEndStateDifference(MiniStep * pMiniStep)
{
	this->lendStateDifferences.push_back(pMiniStep);
}

/**
 *  Clears the vector of endStateDifferences
 */
void Chain::clearEndStateDifferences()
{
	deallocateVector(this->lendStateDifferences);
}

/**
 * Assigns new ordering keys for the ministeps of this chain. The new
 * ordering key of a ministep is its index in the list.
 */
void Chain::resetOrderingKeys()
{
	int key = 0;
	MiniStep * pMiniStep = this->lpFirst;

	while (pMiniStep)
	{
		pMiniStep->orderingKey(key);
		pMiniStep = pMiniStep->pNext();
		key++;
	}
}


/**
 * Updates the "same option" related structures after the insertion of
 * the given ministep.
 */
void Chain::updateSameOptionPointersOnInsert(MiniStep * pNewMiniStep)
{
	MiniStep * pFirstMiniStep =
		this->lfirstMiniStepPerOption[*pNewMiniStep->pOption()];

	if (!pFirstMiniStep ||
		pFirstMiniStep->orderingKey() > pNewMiniStep->orderingKey())
	{
		// The new ministep is the earliest for its option

		pNewMiniStep->pNextWithSameOption(pFirstMiniStep);

		if (pFirstMiniStep)
		{
			pFirstMiniStep->pPreviousWithSameOption(pNewMiniStep);
		}

		this->lfirstMiniStepPerOption[*pNewMiniStep->pOption()] =
			pNewMiniStep;
	}
	else
	{
		MiniStep * pLeftMiniStep = pFirstMiniStep;

		while (pLeftMiniStep->pNextWithSameOption() &&
			pLeftMiniStep->pNextWithSameOption()->orderingKey() <
				pNewMiniStep->orderingKey())
		{
			pLeftMiniStep = pLeftMiniStep->pNextWithSameOption();
		}

		MiniStep * pRightMiniStep = pLeftMiniStep->pNextWithSameOption();

		pLeftMiniStep->pNextWithSameOption(pNewMiniStep);
		pNewMiniStep->pPreviousWithSameOption(pLeftMiniStep);

		if (pRightMiniStep)
		{
			pNewMiniStep->pNextWithSameOption(pRightMiniStep);
			pRightMiniStep->pPreviousWithSameOption(pNewMiniStep);
		}
	}
}


/**
 * Updates the vector of consecutive canceling pairs (CCPs) if the
 * CCP status of the given ministep has changed.
 */
void Chain::updateCCPs(MiniStep * pMiniStep)
{
	if (pMiniStep)
	{
		if (pMiniStep->firstOfConsecutiveCancelingPair() &&
			pMiniStep->consecutiveCancelingPairIndex() == -1)
		{
			pMiniStep->consecutiveCancelingPairIndex(
				this->lccpMiniSteps.size());
			this->lccpMiniSteps.push_back(pMiniStep);
		}
		else if (!pMiniStep->firstOfConsecutiveCancelingPair() &&
			pMiniStep->consecutiveCancelingPairIndex() != -1)
		{
			MiniStep * pLastCCPMiniStep =
				this->lccpMiniSteps[this->lccpMiniSteps.size() - 1];
			int index = pMiniStep->consecutiveCancelingPairIndex();
			this->lccpMiniSteps[index] = pLastCCPMiniStep;
			pLastCCPMiniStep->consecutiveCancelingPairIndex(index);
			this->lccpMiniSteps.pop_back();
			pMiniStep->consecutiveCancelingPairIndex(-1);
		}
	}
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the period whose end observations are connected by this chain.
 */
int Chain::period() const
{
	return this->lperiod;
}


/**
 * Returns the first (dummy) ministep in this chain.
 */
MiniStep * Chain::pFirst() const
{
	return this->lpFirst;
}


/**
 * Returns the last (dummy) ministep in this chain.
 */
MiniStep * Chain::pLast() const
{
	return this->lpLast;
}


/**
 * Returns the initial state of the variables of this chain.
 */
const State * Chain::pInitialState() const
{
	return this->lpInitialState;
}

/**
 * Returns the number of ministeps of this chain excluding the first,
 * but including the last dummy ministep.
 */
int Chain::ministepCount() const
{
	return this->lminiSteps.size();
}


/**
 * Returns the number of diagonal ministeps in this chain.
 */
int Chain::diagonalMinistepCount() const
{
	return this->ldiagonalMiniSteps.size();
}

/**
 * Returns the final (actually the end of the updated set of links)
 * reciprocal rate in the chain.
 */
double Chain::finalReciprocalRate() const
{
	return this->lfinalReciprocalRate;
}

/**
 * Stores the final (actually the end of the updated set of links)
 * reciprocal rate in the chain.
 */
void Chain::finalReciprocalRate(double value)
{
	this->lfinalReciprocalRate = value;
}


/**
 * Returns the sum of reciprocal rates over all non-dummy ministeps.
 */
double Chain::mu() const
{
	return this->lmu;
}


/**
 * Returns the sum of squared reciprocal rates over all non-dummy ministeps.
 */
double Chain::sigma2() const
{
	return this->lsigma2;
}


/**
 * Returns the number of CCPs in this chain.
 */
int Chain::consecutiveCancelingPairCount() const
{
	return this->lccpMiniSteps.size();
}


/**
 * Returns the number of network ministeps with missing observed data
 * at the end of the current period.
 */
int Chain::missingNetworkMiniStepCount() const
{
	return this->lmissingNetworkMiniSteps.size();
}


/**
 * Returns the number of behavior ministeps with missing observed data
 * at the end of the current period.
 */
int Chain::missingBehaviorMiniStepCount() const
{
	return this->lmissingBehaviorMiniSteps.size();
}

/**
 * Returns the vector of ministeps representing differences between the
 * current initial state and the observed data at the start of the
 * current period.
 */
const vector<MiniStep *> & Chain::rInitialStateDifferences() const
{
	return this->linitialStateDifferences;
}
/**
 * Returns the vector of ministeps representing differences between the
 * current end state and the observed data at the end of the
 * current period.
 */
const vector<MiniStep *> & Chain::rEndStateDifferences() const
{
	return this->lendStateDifferences;
}

// ----------------------------------------------------------------------------
// Section: Random draws
// ----------------------------------------------------------------------------

/**
 * Returns a random ministep excluding the first dummy ministep.
 */
MiniStep * Chain::randomMiniStep() const
{
	return this->lminiSteps[nextInt(this->lminiSteps.size())];
}


/**
 * Returns a random diagonal ministep.
 */
MiniStep * Chain::randomDiagonalMiniStep() const
{
	return this->ldiagonalMiniSteps[nextInt(this->ldiagonalMiniSteps.size())];
}


/**
 * Returns a random ministep from the given interval.
 */
MiniStep * Chain::randomMiniStep(MiniStep * pFirstMiniStep,
	MiniStep * pLastMiniStep) const
{
	int length = this->intervalLength(pFirstMiniStep, pLastMiniStep);
	int index = nextInt(length);
	MiniStep * pMiniStep = pFirstMiniStep;

	while (index > 0)
	{
		pMiniStep = pMiniStep->pNext();
		index--;
	}

	return pMiniStep;
}


/**
 * Returns the first mini step of a random consecutive canceling pair.
 */
MiniStep * Chain::randomConsecutiveCancelingPair() const
{
	return this->lccpMiniSteps[nextInt(this->lccpMiniSteps.size())];
}


/**
 * Returns a random network change with missing observed data.
 */
MiniStep * Chain::randomMissingNetworkMiniStep() const
{
	return this->lmissingNetworkMiniSteps[
		nextInt(this->lmissingNetworkMiniSteps.size())];
}


/**
 * Returns a random behavior change with missing observed data.
 */
MiniStep * Chain::randomMissingBehaviorMiniStep() const
{
	return this->lmissingBehaviorMiniSteps[
		nextInt(this->lmissingBehaviorMiniSteps.size())];
}


// ----------------------------------------------------------------------------
// Section: Intervals
// ----------------------------------------------------------------------------

/**
 * Returns the length of the given interval of ministeps.
 */
int Chain::intervalLength(const MiniStep * pFirstMiniStep,
	const MiniStep * pLastMiniStep) const
{
	int length = 1;
	const MiniStep * pMiniStep = pFirstMiniStep;

	while (pMiniStep != pLastMiniStep)
	{
		pMiniStep = pMiniStep->pNext();
		length++;
	}

	return length;
}


// ----------------------------------------------------------------------------
// Section: Same option related
// ----------------------------------------------------------------------------

/**
 * Returns the first ministep of the given option or 0, if there is no such a
 * ministep.
 */
MiniStep * Chain::firstMiniStepForOption(const Option & rOption) const
{
	MiniStep * pMiniStep = 0;
	map<const Option, MiniStep *>::const_iterator iter =
		this->lfirstMiniStepPerOption.find(rOption);

	if (iter != this->lfirstMiniStepPerOption.end())
	{
		pMiniStep = iter->second;
	}

	return pMiniStep;
}


/**
 * Returns the first ministep of the given option in the subchain starting
 * with the given ministep, or 0, if there is no such a ministep.
 */
MiniStep * Chain::nextMiniStepForOption(const Option & rOption,
	const MiniStep * pFirstMiniStep) const
{
	MiniStep * pMiniStep = this->firstMiniStepForOption(rOption);

	while (pMiniStep &&
		pMiniStep->orderingKey() < pFirstMiniStep->orderingKey())
	{
		pMiniStep = pMiniStep->pNextWithSameOption();
	}

	return pMiniStep;
}
// ----------------------------------------------------------------------------
// Section: Same dyad related
// ----------------------------------------------------------------------------

/**
 * Returns the first ministep of the ego, alter, actor set1, actor set 2,
 * or pLast if there is no such a ministep.
 */
MiniStep * Chain::pFirstMiniStepForLink(const MiniStep * pLinkMiniStep) const
{
	int ego = pLinkMiniStep->ego();
	const NetworkChange * pLinkNetworkChange =
		dynamic_cast<const NetworkChange *>(pLinkMiniStep);
	int alter = pLinkNetworkChange->alter();
	const ActorSet * pActorSet1 =
		this->lpData->pNetworkData(pLinkMiniStep->variableName())->pSenders();
	const ActorSet * pActorSet2 =
		this->lpData->pNetworkData(pLinkMiniStep->variableName())->pReceivers();
	bool done = false;
	MiniStep * pMiniStep = this->lpFirst->pNext();
	while (pMiniStep != this->lpLast && !done)
	{
		if (pMiniStep->networkMiniStep())
		{
			const NetworkChange * pNetworkChange =
				dynamic_cast<const NetworkChange *>(pMiniStep);
			if (pMiniStep->ego() == ego &&
				pNetworkChange->alter() == alter)
			{
				NetworkLongitudinalData * pNetwork =
					this->lpData->pNetworkData(pMiniStep->variableName());
				if (pNetwork->pSenders() == pActorSet1 &&
					pNetwork->pReceivers() == pActorSet2)
				{
					done = true;
					break;
				}
			}
		}
		pMiniStep = pMiniStep->pNext();
	}
	if (pMiniStep != this->lpLast)
	{
		Rf_PrintValue(getMiniStepDF(*pMiniStep));
	}
	else
		Rprintf("last\n");
	return pMiniStep;
}


/**
 * Returns the last ministep of the given option in the subchain starting
 * with the given ministep, or pFirst, if there is no such a ministep.
 */
MiniStep * Chain::pLastMiniStepForLink(const MiniStep * pLinkMiniStep) const
{
	int ego = pLinkMiniStep->ego();
	const NetworkChange * pLinkNetworkChange =
		dynamic_cast<const NetworkChange *>(pLinkMiniStep);
	int alter = pLinkNetworkChange->alter();
	const ActorSet * pActorSet1 =
		this->lpData->pNetworkData(pLinkMiniStep->variableName())->pSenders();
	const ActorSet * pActorSet2 =
		this->lpData->pNetworkData(pLinkMiniStep->variableName())->pReceivers();
	MiniStep * pLastMiniStep = this->lpFirst;
	MiniStep * pMiniStep = this->lpFirst->pNext();

	while (pMiniStep != this->lpLast)
	{
		if (pMiniStep->networkMiniStep())
		{
			const NetworkChange * pNetworkChange =
				dynamic_cast<const NetworkChange *>(pMiniStep);
			if (pMiniStep->ego() == ego &&
				pNetworkChange->alter() == alter)
			{
				NetworkLongitudinalData * pNetwork =
					this->lpData->pNetworkData(pMiniStep->variableName());
				if (pNetwork->pSenders() == pActorSet1 &&
					pNetwork->pReceivers() == pActorSet2)
				{
					pLastMiniStep = pMiniStep;
				}
			}
		}
		pMiniStep = pMiniStep->pNext();
	}
	return pLastMiniStep;
}
// ----------------------------------------------------------------------------
// Section: copy a chain (too difficult (for Ruth!) to do a copy constructor!
// ----------------------------------------------------------------------------

/**
 * Returns a copy of the given chain
 */
Chain * Chain::copyChain() const
{
	Chain * pChain = new Chain(this->lpData);

	pChain->lperiod = this->lperiod;

	MiniStep *pMiniStep = this->lpFirst->pNext();
	while (pMiniStep != this->lpLast)
//	for (unsigned i = 1; i < this->lminiSteps.size(); i++)
	{

		MiniStep *pCopyMiniStep = pMiniStep->createCopyMiniStep();

		// need the reciprocal rate to get the mu and sigma2 correct
		pCopyMiniStep->reciprocalRate(pMiniStep->reciprocalRate());

		pChain->insertBefore(pCopyMiniStep, pChain->lpLast);
		pMiniStep = pMiniStep->pNext();
	}

	pChain->lmu = this->lmu;
	pChain->lsigma2 = this->lsigma2;
//	Rprintf("%x\n", pChain);
	//	SEXP ch2 = getChainDF(*this);
				//		Rf_PrintValue(ch2);
				//	SEXP ch1 = getChainDF(*pChain);
				//		Rf_PrintValue(ch1);
	for(unsigned i=0; i < this->linitialStateDifferences.size(); i++)
	{
		pChain->linitialStateDifferences.push_back(
			this->linitialStateDifferences[i]->createCopyMiniStep());
	}
	for(unsigned i=0; i < this->lendStateDifferences.size(); i++)
	{
		pChain->lendStateDifferences.push_back(
			this->lendStateDifferences[i]->createCopyMiniStep());
	}
// 	for(int i=0; i< this->linitialStateDifferences.size(); i++)
// 	{
// 		Rf_PrintValue(getMiniStepDF(*(this->linitialStateDifferences[i])));
// 		Rf_PrintValue(getMiniStepDF(*(pChain->linitialStateDifferences[i])));
// 	}
	return pChain;
	}

//void Chain::dumpChain() const
// see RSienaTest


}

