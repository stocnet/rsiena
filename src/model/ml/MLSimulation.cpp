/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MLSimulation.cpp
 *
 * Description: This file contains the implementation of the class
 * MLSimulation.
 * The Metropolis-Hastings steps are extensively documented
 * in Siena_algorithms.pdf
 *****************************************************************************/
#include <stdexcept>
#include <string>
#include <cmath>
#include <R_ext/Error.h>
#include "MLSimulation.h"

#include <Rinternals.h>
#include "utils/Random.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "network/TieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/Model.h"
#include "model/State.h"
#include "model/SimulationActorSet.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/BehaviorChange.h"
#include "model/tables/Cache.h"
#include <Rinternals.h>

using namespace std;

namespace siena
{
SEXP getMiniStepDF(const MiniStep& miniStep);
SEXP getChainDF(const Chain& chain, bool sort=true);

// ----------------------------------------------------------------------------
// Section: Constructors and destructors
// ----------------------------------------------------------------------------

/**
 * Creates a MLsimulation object for the given observed data and an
 * actor-based model for that data.
 */
MLSimulation::MLSimulation(Data * pData, Model * pModel) :
	EpochSimulation(pData, pModel)
{
	for (int i = 0; i < NBRTYPES; i++)
	{
//		this->lrejections[i] = 0;
//		this->lacceptances[i] = 0;
		this->laborts[i] = 0;
	}
	this->lthisPermutationLength = 0;
}


/**
 * Deallocates this simulation object.
 */
MLSimulation::~MLSimulation()
{
	deallocateVector(this->linitialMissingOptions);
}

// ----------------------------------------------------------------------------
// Section: Model simulations
// ----------------------------------------------------------------------------

/**
 * Initializes the simulation as of the beginning of the specified
 * period.
 */
void MLSimulation::initialize(int period)
{
	EpochSimulation::initialize(period);

	deallocateVector(this->linitialMissingOptions);

	for (unsigned k = 0;
		k < this->pData()->rDependentVariableData().size();
		k++)
	{
		const NetworkLongitudinalData * pNetworkData =
			dynamic_cast<const NetworkLongitudinalData *>(
				this->pData()->rDependentVariableData()[k]);
		const BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<const BehaviorLongitudinalData *>(
				this->pData()->rDependentVariableData()[k]);

		if (pNetworkData)
		{
			for (TieIterator iter =
					pNetworkData->pMissingTieNetwork(period)->ties();
				iter.valid();
				iter.next())
			{
				this->linitialMissingOptions.push_back(
					new Option(pNetworkData->id(), iter.ego(), iter.alter()));
			}
		}
		else if (pBehaviorData)
		{
			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				if (pBehaviorData->missing(period, i))
				{
					this->linitialMissingOptions.push_back(
						new Option(pBehaviorData->id(), i));
				}
			}
		}
	}
}

/**
 * Initializes the initial state from the observed data for the specified
 * period. Not to be done too often. Once a chain exists it needs its
 * initial state till it is destroyed without further initialization.
 */
void MLSimulation::initializeInitialState(int period)
{
	// Create the initial state. The observed values in the Data object
	// should be copied if there are any missings.
	bool copyValues = !this->linitialMissingOptions.empty();
	this->pChain()->setupInitialState(copyValues);
}


/**
 * Generates a random chain connecting the start and end observations of the
 * given data object for the given period. The chain is simple in the sense
 * that no two ministeps cancel each other out.
 */
void MLSimulation::connect(int period)
{
	this->initialize(period);
	this->initializeInitialState(period);
	this->pChain()->connect(period, this);
	this->updateProbabilities(this->pChain(), this->pChain()->pFirst()->pNext(),
			this->pChain()->pLast()->pPrevious());
}

/**
 * Updates the probabilities for a range of ministeps of the given chain
 * (including the end-points of the range).
 */
void MLSimulation::updateProbabilities(Chain * pChain,
	MiniStep * pFirstMiniStep,
	MiniStep * pLastMiniStep)
{
	// Initialize the variables as of the beginning of the period
	resetVariables();

    // Apply the ministeps before the first ministep of the required range
    // to derive the correct state before that ministep.

    this->executeMiniSteps(pChain->pFirst()->pNext(), pFirstMiniStep);

    bool done = false;
    MiniStep * pMiniStep = pFirstMiniStep;

	// set up array to store counts of structurally active ministeps by variable
	int *counts = new int[this->lvariables.size()];
	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
		counts[i] = 0;
	}
	while (!done)
    {
    	DependentVariable * pVariable =
			this->lvariables[pMiniStep->variableId()];
		this->calculateRates();
		double rate = pVariable->rate(pMiniStep->ego());
		double probability = pVariable->probability(pMiniStep);
		double reciprocalTotalRate = 1 / this->grandTotalRate();

		if (!pVariable->structural(pMiniStep))
		{
			counts[pMiniStep->variableId()] ++;
		}
		//Rprintf("ValidMinistep %d\n", pVariable->validMiniStep(pMiniStep));
		//Rf_PrintValue(getMiniStepDF(*pMiniStep));
		pMiniStep->reciprocalRate(reciprocalTotalRate);
		pMiniStep->logOptionSetProbability(log(rate * reciprocalTotalRate));
		pMiniStep->logChoiceProbability(log(probability));
		pMiniStep->makeChange(pVariable);

    	if (pMiniStep == pLastMiniStep)
    	{
    		done = true;
    	}
    	else
    	{
    		pMiniStep = pMiniStep->pNext();
    	}
    }
	// calculate rate scores for the dependent variables
	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
		this->lvariables[i]->calculateMaximumLikelihoodRateScores(counts[i]);
	}

	// if necessary, calculate rate derivatives for the dependent variables
	if(this->pModel()->needDerivatives())
	{
		for (unsigned i = 0; i < this->lvariables.size(); i++)
		{
			this->lvariables[i]->
				calculateMaximumLikelihoodRateDerivatives(counts[i]);
		}

	}
	delete [] counts;
	// finally (sometimes) we need the next reciprocal rate.
	this->calculateRates();
	pChain->finalReciprocalRate(1 / this->grandTotalRate());

}

/**
 *  Does a few steps to increase the length of a newly constructed
 *  minimal chain.
 */

void MLSimulation::preburnin()
{
	int rejectCount = 0;
	bool accept;
	while (rejectCount < 50)
	{
		accept = this->insertDiagonalMiniStep();
		if (!accept)
		{
			rejectCount ++;
		}
	}
	rejectCount = 0;
	while (rejectCount < 50)
	{
		accept = this->insertPermute(1);
		if (!accept)
		{
			rejectCount ++;
		}
	}
}


/**
 *  Merges the probabilities, and does the required number of steps
 *  in the current chain
 *
 */

void MLSimulation::runEpoch(int period)
{
	// get the variables from the observed values
	this->initialize(period);
	// create an initial state object
	this->initializeInitialState(period);
	// copy over the initial state differences from the stored vector
	this->pChain()->recreateInitialState();

	this->setUpProbabilityArray();

	// recreate the probabilities on the chain using the current parameters
	this->updateProbabilities(this->pChain(),
			this->pChain()->pFirst()->pNext(),
			this->pChain()->pLast()->pPrevious());

	int numSteps = this->pModel()->numberMLSteps() ;

	for (int i = 0; i < numSteps; i++)
	{
		this->MLStep();
	}
}

/**
 *  sets up probability array
 *
 */
void MLSimulation::setUpProbabilityArray()
{
	this->lprobabilityArray[0] = this->pModel()->insertDiagonalProbability();
	this->lprobabilityArray[1] = this->pModel()->cancelDiagonalProbability();
	this->lprobabilityArray[2] = this->pModel()->permuteProbability();
	this->lprobabilityArray[3] = this->pModel()->insertPermuteProbability();
	this->lprobabilityArray[4] = this->pModel()->deletePermuteProbability();
	this->lprobabilityArray[5] =
		this->pModel()->insertRandomMissingProbability();
	this->lprobabilityArray[6] =
		this->pModel()->deleteRandomMissingProbability();
	this->lprobabilityArray[7] = 1 - this->lprobabilityArray[0] -
		this->lprobabilityArray[1] - this->lprobabilityArray[2] -
		this->lprobabilityArray[3] - this->lprobabilityArray[4] -
		this->lprobabilityArray[5] - this->lprobabilityArray[6];

	for (int i = 0; i < NBRTYPES; i++)
	{
//		this->lrejections[i] = 0;
//		this->lacceptances[i] = 0;
		this->laborts[i] = 0;
	}
}

/**
 *  Does one step in the current chain
 *
 */
void MLSimulation::MLStep()
{
	int stepType = nextIntWithProbabilities(NBRTYPES, this->lprobabilityArray);
	int c0 = this->lcurrentPermutationLength;
//	Rf_PrintValue(getChainDF(*this->pChain()));
	bool accept = false;
	this->lproposalProbability = R_NaN;

	switch (stepType)
	{
	case 0:
		accept = this->insertDiagonalMiniStep();
		break;
	case 1:
		accept = this->cancelDiagonalMiniStep();
		break;
	case 2:
		accept = this->permute(c0);
		this->updateCurrentPermutationLength(accept);
		break;
	case 3:
		accept = this->insertPermute(c0);
		this->updateCurrentPermutationLength(accept);
		break;
	case 4:
		accept = this->deletePermute(c0);
		this->updateCurrentPermutationLength(accept);
		break;
	case 5:
		accept = this->insertMissing();
		break;
	case 6:
		accept = this->deleteMissing();
		break;
	case 7:
		accept = this->move();
		break;
	}
// The following is done within the classes insertDiagonalMiniStep etc.
// as part of recordOutcome
//	if (accept)
//	{
//		this->lacceptances[stepType]++;
//	}
//	else if (!R_IsNaN(this->lproposalProbability))
//	{
//		this->lrejections[stepType]++;
//	}
//	else
//	{
//		this->laborted[stepType]++;
//	}
//	Rprintf("%d %d %f\n",stepType, accept, this->lproposalProbability);
}

/**
 * records the outcome of the current MH step
 */
void MLSimulation::recordOutcome(const MiniStep & miniStep, bool accept,
	int stepType, bool misdat)
{
	DependentVariable * pVariable =
		this->lvariables[miniStep.variableId()];

	if (stepType == PERMUTE)
	{
		pVariable =	this->lvariables[0];
	}
	if (misdat) // changes stepType INSPERM or DELPERM to INSMISDAT, DELMISDAT.
	{
		stepType = stepType + 4;
	}
	if (accept)
	{
		pVariable->incrementAcceptances(stepType);
	}
	else if (!R_IsNaN(this->lproposalProbability))
	{
		pVariable->incrementRejections(stepType);
	}
	else
	{
//		pVariable->incrementAborts(stepType);
		this->incrementAborts(stepType);
	}
}

/**
 * sets the current state to immediately before the given mini step.
 */
void MLSimulation::setStateBefore(MiniStep * pMiniStep)
{
	this->resetVariables();
	this->executeMiniSteps(this->pChain()->pFirst()->pNext(), pMiniStep);
}

/**
 * Executes the given subsequence of ministeps excluding the last
 * ministep.
 */
void MLSimulation::executeMiniSteps(MiniStep * pFirstMiniStep,
	MiniStep * pLastMiniStep)
{
	MiniStep * pMiniStep = pFirstMiniStep;

	while (pMiniStep != pLastMiniStep)
	{
		DependentVariable * pVariable =
			this->lvariables[pMiniStep->variableId()];
		pMiniStep->makeChange(pVariable);

		pMiniStep = pMiniStep->pNext();
	}
}

/**
 * Resets the variables when it is necessary to read the initial state.
 */
void MLSimulation::resetVariables()
{
	// Initialize each dependent variable

	for (unsigned i = 0; i < this->rVariables().size(); i++)
	{
		DependentVariable * pVariable = this->rVariables()[i];
		pVariable->initialize(this->period());

		// The values of missings have to be read from the current
		// initial state (y_init in the specification).
		if (!this->linitialMissingOptions.empty())
		{
			if (pVariable->networkVariable())
			{
				const Network * pInitialNetwork =
					this->pChain()->pInitialState()->pNetwork(
						pVariable->name());
				NetworkVariable * pNetworkVariable =
					dynamic_cast<NetworkVariable *>(pVariable);
				NetworkLongitudinalData * pNetworkData =
					dynamic_cast<NetworkLongitudinalData *>(
						pNetworkVariable->pData());
				const Network * pMissings =
					pNetworkData->pMissingTieNetwork(this->period());

				for (TieIterator iter = pMissings->ties();
					iter.valid();
					iter.next())
				{
					pNetworkVariable->pNetwork()->setTieValue(iter.ego(),
						iter.alter(),
						pInitialNetwork->tieValue(iter.ego(),
							iter.alter()));
				}
			}
			else if (pVariable->behaviorVariable())
			{
				const int * initialValues =
					this->pChain()->pInitialState()->behaviorValues(
						pVariable->name());
				BehaviorVariable * pBehaviorVariable =
					dynamic_cast<BehaviorVariable *>(pVariable);
				BehaviorLongitudinalData * pBehaviorData =
					dynamic_cast<BehaviorLongitudinalData *>(
						pBehaviorVariable->pData());

				for (int actor = 0; actor < pBehaviorData->n(); actor++)
				{
					if (pBehaviorData->missing(this->period(), actor))
					{
						pBehaviorVariable->value(actor, initialValues[actor]);
					}
				}
			}
		}
	}
}

// ----------------------------------------------------------------------------
// Section: Metropolis-Hastings steps
// ----------------------------------------------------------------------------

/**
 * Implements the MH_InsDiag Metropolis-Hastings step.
 */
bool MLSimulation::insertDiagonalMiniStep()
{
	bool accept = false;
	MiniStep * pMiniStep = this->pChain()->randomMiniStep();
	this->setStateBefore(pMiniStep);
	this->calculateRates();
	DependentVariable * pVariable = this->chooseVariable();
	int i = this->chooseActor(pVariable);
	BehaviorVariable * pBehaviorVariable =
		dynamic_cast<BehaviorVariable *>(pVariable);
	NetworkVariable * pNetworkVariable =
		dynamic_cast<NetworkVariable *>(pVariable);
	if (!pVariable->pActorSet()->active(i) ||
		(pBehaviorVariable && pBehaviorVariable->structural(i)))
	{
		return false;
	}

	MiniStep * pNewMiniStep = 0;

	if (pBehaviorVariable)
	{
		pNewMiniStep =
			new BehaviorChange(
				dynamic_cast<BehaviorLongitudinalData *>(pVariable->pData()),
				i,
				0);
	}
	else
	{
		if (pNetworkVariable->oneModeNetwork())
		{
			pNewMiniStep =
				new NetworkChange(
					dynamic_cast<NetworkLongitudinalData *>(pVariable->pData()),
					i,
					i, true);
		}
		else
		{
			pNewMiniStep =
				new NetworkChange(
					dynamic_cast<NetworkLongitudinalData *>(pVariable->pData()),
					i,
					pVariable->m(), true);

		}
	}

	double rr = 1 / this->grandTotalRate();
	pNewMiniStep->reciprocalRate(rr);
	pNewMiniStep->logOptionSetProbability(log(pVariable->rate(i) * rr));

	double explcpr = pVariable->probability(pNewMiniStep);
	pNewMiniStep->logChoiceProbability(log(explcpr));

	double kappaFactor;

	if (this->simpleRates())
	{
		kappaFactor = 1 / (rr * (this->pChain()->ministepCount()));
	}
	else
	{
		double sigma2 = this->pChain()->sigma2();
		double mu = this->pChain()->mu();

		kappaFactor = sqrt(sigma2 / (sigma2 + rr * rr)) *
			exp((1 - mu) * (1 - mu) / (2 * sigma2) -
				(1 - mu - rr) * (1 - mu - rr) / (2 * (sigma2 + rr * rr)));
	}

	this->lproposalProbability =
		kappaFactor *
		explcpr *
		this->pChain()->ministepCount() *
		this->pModel()->cancelDiagonalProbability() /
		((this->pChain()->diagonalMinistepCount() + 1) *
			this->pModel()->insertDiagonalProbability());

//	if (this->lproposalProbability > 1)
//	{
//		this->lproposalProbability = 1;
//	}

	accept = nextDouble() < this->lproposalProbability;
	if (accept)
	{
		this->pChain()->insertBefore(pNewMiniStep, pMiniStep);
	}

	this->recordOutcome(*pNewMiniStep, accept, INSDIAG, false);
	if (!accept)
	{
		delete pNewMiniStep;
	}
// don't I need to delete the ministep if rejected?
	return accept;
}

/**
 * Implements the MH_CancDiag Metropolis-Hastings step.
 */
bool MLSimulation::cancelDiagonalMiniStep()
{
	bool accept = false;

	if (this->pChain()->diagonalMinistepCount() == 0)
	{
		return accept;
	}

	MiniStep * pMiniStep = this->pChain()->randomDiagonalMiniStep();
	double rr = pMiniStep->reciprocalRate();

	double kappaFactor;

	if (this->simpleRates())
	{
		kappaFactor = rr * (this->pChain()->ministepCount() - 1);
	}
	else
	{
		double sigma2 = this->pChain()->sigma2();
		double mu = this->pChain()->mu();

		kappaFactor = sqrt(sigma2 / (sigma2 + rr * rr)) *
			exp((1 - mu) * (1 - mu) / (2 * sigma2) -
				(1 - mu + rr) * (1 - mu + rr) / (2 * (sigma2 - rr * rr)));
	}

	this->lproposalProbability =
		kappaFactor * exp(-pMiniStep->logChoiceProbability()) *
			this->pChain()->diagonalMinistepCount() *
		this->pModel()->insertDiagonalProbability() /
		((this->pChain()->ministepCount() - 1) *
			this->pModel()->cancelDiagonalProbability());

	if (this->lproposalProbability > 1)
	{
		this->lproposalProbability = 1;
	}

	accept = nextDouble() < this->lproposalProbability;

	this->recordOutcome(*pMiniStep, accept, CANCDIAG, false);

	if (accept)
	{
		this->pChain()->remove(pMiniStep);
		delete pMiniStep;
	}

	return accept;
}

/**
 * Implements the MH_move Metropolis-Hastings step.
 */
bool MLSimulation::move()
{
	//// not implemented for non simple rates
	if (!this->simpleRates())
	{
		return false;
	}

	if (this->pChain()->ministepCount() <= 2)
	{
		return false;
	}

	bool accept = false;

	MiniStep * pMiniStepA = this->pChain()->randomMiniStep();

	while ((pMiniStepA == this->pChain()->pLast()) ||
		(pMiniStepA->diagonal()))
	{
		pMiniStepA = this->pChain()->randomMiniStep();
	}

	// We will move A to before a ministep in the interval {pFirst, ..., pLast}
	MiniStep * pFirst = this->pChain()->pFirst()->pNext();
	MiniStep * pLast = this->pChain()->pLast();

	MiniStep * pMiniStepD = pMiniStepA->pNextWithSameOption();

	if (pMiniStepD)
	{
		if (pMiniStepD->pPrevious()==pMiniStepA)
		{
			return false;
		}
		pLast = pMiniStepD->pPrevious();
	}

	MiniStep * pMiniStepE = pMiniStepA->pPreviousWithSameOption();

	if (pMiniStepE)
	{
		if (pMiniStepE->pNext()==pMiniStepA)
		{
			return false;
		}
		pFirst = pMiniStepE->pNext();
		pFirst = pFirst->pNext();
	}

	int intLength = this->pChain()->intervalLength(pFirst, pLast)-1;

	//// Position of ministepA

	int piCount = 0;
	MiniStep * pMiniStepC = pFirst;
	while (pMiniStepC != pMiniStepA)
	{
		piCount++;
		pMiniStepC = pMiniStepC->pNext();
	}
	int positionA = piCount;

	////CHOOSING MINISTEPB

	double * newReciprocalRate = new double[intLength];
	double * newOptionSetProbability = new double[intLength];
	double * newChoiceProbability = new double[intLength];

	double * reciprocalRateA = 0;
	double * optionSetProbabilityA = 0;
	double * choiceProbabilityA = 0;

	if (this->pModel()->localML())
	{
		reciprocalRateA = new double[intLength];
		optionSetProbabilityA = new double[intLength];
		choiceProbabilityA = new double[intLength];
	}

	bool * updateProb = new bool[intLength];
	double * pi = new double[intLength];

	pi[0] = 1;
	piCount = 0;

	// SETTING PI AND UPDATE
	for (int k=0;k < intLength;k++)
	{
		updateProb[k] = false;
	}
	for (int k=1;k < intLength;k++)
	{
		pi[k] = 0;
	}

	MiniStep * pMiniStepAReverse = pMiniStepA->createReverseMiniStep();
	pMiniStepC = pFirst;

	// CREATING VARIABLES, EGOS AND ALTERS FOR MiniSteps A AND C
	DependentVariable * pVariableA = this->lvariables[pMiniStepA->variableId()];
	NetworkVariable * pNetworkVariableA =
		dynamic_cast<NetworkVariable *>(pVariableA);
	DependentVariable * pVariableC = this->lvariables[pMiniStepC->variableId()];
	NetworkVariable * pNetworkVariableC =
		dynamic_cast<NetworkVariable *>(pVariableC);
	int egoA = pMiniStepA->ego();
	int alterA = 0;
	if (pMiniStepA->networkMiniStep())
	{
		alterA = dynamic_cast<NetworkChange *>(pMiniStepA)->alter();
	}
	int egoC = pMiniStepC->ego();
	int alterC = 0;
	if (pMiniStepC->networkMiniStep())
	{
		alterC = dynamic_cast<NetworkChange *>(pMiniStepC)->alter();
	}

	this->setStateBefore(pFirst);

	this->calculateRates();
	double rrA = 1 / this->grandTotalRate();
	double osprA = pVariableA->rate(egoA) * rrA;
	double cprA = pVariableA->probability(pMiniStepA);
	double rHere = osprA * cprA;

	if (this->pModel()->localML())
	{
		reciprocalRateA[0] = rrA;
		optionSetProbabilityA[0] = log(osprA);
		choiceProbabilityA[0] = log(cprA);
	}

	bool Adone = false;

	// Calculating values of pi

	while (pMiniStepC != pLast)
	{
		if (pMiniStepC == pMiniStepA)
		{
			Adone = true;
			pMiniStepC->makeChange(pVariableC);
		} else {
			piCount++;
			pi[piCount] = pi[piCount-1];

			if (pVariableC->pActorSet() == pVariableA->pActorSet())
			{
				bool N2 = true;
				if (this->pModel()->localML())
				{
					N2 = this->neighbourhoodChange(pMiniStepA,pMiniStepC,
						pVariableA,pNetworkVariableA,egoA,alterA);
				} else {
					N2 = this->smallNeighbourhoodChange(pMiniStepA,pMiniStepC,
														pVariableA,pNetworkVariableA,egoA,alterA);
				}
				if (N2)
				{
					double rr = 0;
					double ospr = 0;
					double cpr =  0;
					double cCurrent = 0;

					if (Adone)
					{
						pMiniStepAReverse->makeChange(pVariableA);
						this->calculateRates();
						rr = 1 / this->grandTotalRate();
						ospr = pVariableC->rate(egoC) * rr;
						cpr = pVariableC->probability(pMiniStepC);
						pMiniStepA->makeChange(pVariableA);
						cCurrent = pMiniStepC->logChoiceProbability() +
							pMiniStepC->logOptionSetProbability();
						pi[piCount] *= ((cpr*ospr)/exp(cCurrent));

					} else {
						pMiniStepA->makeChange(pVariableA);
						this->calculateRates();
						rr = 1 / this->grandTotalRate();
						ospr = pVariableC->rate(egoC) * rr;
						cpr = pVariableC->probability(pMiniStepC);
						pMiniStepAReverse->makeChange(pVariableA);
						cCurrent = pMiniStepC->logChoiceProbability() +
							pMiniStepC->logOptionSetProbability();
						pi[piCount] *= (exp(cCurrent)/(cpr*ospr));
					}
					newReciprocalRate[piCount] = rr;
					newOptionSetProbability[piCount] = log(ospr);
					newChoiceProbability[piCount] = log(cpr);
					updateProb[piCount] = true;
				}

				pMiniStepC->makeChange(pVariableC);

				if (Adone) {pMiniStepAReverse->makeChange(pVariableA);}

				N2 = true;
				if (this->pModel()->localML())
				{
					N2 = this->neighbourhoodChange(pMiniStepC,pMiniStepA,
						pVariableC,
						pNetworkVariableC,
						egoC,alterC);
				} else {
					N2 = this->smallNeighbourhoodChange(pMiniStepC,pMiniStepA,
														pVariableC,
														pNetworkVariableC,
														egoC,alterC);
				}

				if (N2)
				{
					double rPrev = rHere;
					this->calculateRates();
					rrA = 1 / this->grandTotalRate();
					osprA = pVariableA->rate(egoA) * rrA;
					cprA = pVariableA->probability(pMiniStepA);
					rHere = osprA * cprA;
					pi[piCount] *= (rHere/rPrev);
				}
				if (Adone) {pMiniStepA->makeChange(pVariableA);}

			} else {//if actor set don't match
				pMiniStepC->makeChange(pVariableC);
			}

			if (this->pModel()->localML())
			{
				reciprocalRateA[piCount] = rrA;
				optionSetProbabilityA[piCount] = log(osprA);
				choiceProbabilityA[piCount] = log(cprA);
			}
		}

		pMiniStepC = pMiniStepC->pNext();
		if (pMiniStepC != this->pChain()->pLast())
		{
			pVariableC = this->lvariables[pMiniStepC->variableId()];
			pNetworkVariableC =
				dynamic_cast<NetworkVariable *>(pVariableC);
			egoC = pMiniStepC->ego();
			alterC = 0;
			if (pMiniStepC->networkMiniStep()) {alterC = dynamic_cast<NetworkChange *>(pMiniStepC)->alter();}
		}
	}

	delete pMiniStepAReverse;

	// Making pi sum to 1
	double piSum = 0;
	for (int k=0;k<intLength;k++)
	{
		piSum += pi[k];
	}
	for (int k=0;k<intLength;k++)
	{
		pi[k] /= piSum;
	}

	int miniStepChoice = nextIntWithProbabilities(intLength, pi);

	// Will give no change
		if (miniStepChoice == positionA)
		{
			return false;
		}

	double lcprA = 0;
	double losprA = 0;

	if (this->pModel()->localML())
	{
		accept = true;
		lcprA = choiceProbabilityA[miniStepChoice];
		losprA = optionSetProbabilityA[miniStepChoice];
		rrA = reciprocalRateA[miniStepChoice];
		this->lproposalProbability = 1;

	} else {
		MiniStep * pMiniStepB = pFirst;

		for (int k=0;k<miniStepChoice;k++)
		{
			pMiniStepB = pMiniStepB->pNext();
		}
		if (!(miniStepChoice < positionA))
		{
			pMiniStepB = pMiniStepB->pNext();
		}
	// A proposed to be inserted before B
		double logProbChange = -pMiniStepA->logChoiceProbability();

		MiniStep * pFirstA = 0;
		MiniStep * pSecondA = 0;

		int count = 0;
		if (miniStepChoice < positionA)
		{
			pFirstA = pMiniStepB;
			pSecondA = pMiniStepA;
			count = miniStepChoice;
		} else {
			pFirstA = pMiniStepA;
			pSecondA = pMiniStepB;
			count = positionA;
		}

		double rr = 0;
		double lcpr = 0;
		double ospr = 0;
		double cpr = 0;
		this->setStateBefore(pFirstA);
		if (miniStepChoice < positionA)
		{
			this->calculateRates();
			rrA = 1 / this->grandTotalRate();
			osprA = pVariableA->rate(egoA) * rrA;
			losprA = log(osprA);
			lcprA = log(pVariableA->probability(pMiniStepA));
			logProbChange += lcprA;
			pMiniStepA->makeChange(pVariableA);
			pMiniStepC = pFirstA;
		} else {
			pMiniStepC = pFirstA->pNext();
		}
		pVariableC = this->lvariables[pMiniStepC->variableId()];

		// Calculate probabilities between pFirstA and pSecondA
		while (pMiniStepC != pSecondA)
		{
			count++;
			if (pVariableC->pActorSet() == pVariableA->pActorSet())
			{
				if (updateProb[count])
				{
					lcpr = newChoiceProbability[count];

				} else {
					this->calculateRates();
					rr = 1 / this->grandTotalRate();
					cpr = pVariableC->probability(pMiniStepC);
					ospr = pVariableC->rate(egoC) * rr;
					lcpr = log(cpr);
					newReciprocalRate[count] = rr;
					newOptionSetProbability[count] = log(ospr);
					newChoiceProbability[count] = lcpr;
					updateProb[count] = true;
				}
				logProbChange += lcpr;
				logProbChange -= pMiniStepC->logChoiceProbability();
			}

			pMiniStepC->makeChange(pVariableC);
			pMiniStepC = pMiniStepC->pNext();
			if (pMiniStepC != this->pChain()->pLast())
			{
				pVariableC = this->lvariables[pMiniStepC->variableId()];
			}
		}

		if (!(miniStepChoice < positionA))
		{
			this->calculateRates();
			rrA = 1 / this->grandTotalRate();
			osprA = pVariableA->rate(egoA) * rrA;
			losprA = log(osprA);
			lcprA = log(pVariableA->probability(pMiniStepA));
			logProbChange += lcprA;
		}

		this->lproposalProbability = exp(logProbChange)*pi[positionA]/pi[miniStepChoice];
		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}
	}

	accept = nextDouble() < this->lproposalProbability;

	if (accept)
	{
	// Move miniStepA, and update the necessary probabilities
		pMiniStepC = pFirst;

	//Otherwise pMiniStepC will be removed below.
		if (pFirst == pMiniStepA)
		{
			pMiniStepC = pMiniStepC->pNext();
		}
		this->pChain()->remove(pMiniStepA);

		bool startUpdating = false;
		bool stopUpdating = false;
		piCount = 0;

		while(pMiniStepC != pLast)
		{
			if (piCount == miniStepChoice)
			{
				this->pChain()->insertBefore(pMiniStepA, pMiniStepC);
				stopUpdating = startUpdating;
				startUpdating = true;
			}
			if (piCount == positionA)
			{
				stopUpdating = startUpdating;
				startUpdating = true;
			}
			if (stopUpdating)
			{
				break;
			}
			piCount++;

			if (startUpdating)
			{
				if (updateProb[piCount])
				{
					pMiniStepC->logChoiceProbability(newChoiceProbability[piCount]);
					pMiniStepC->reciprocalRate(newReciprocalRate[piCount]);
					pMiniStepC->logOptionSetProbability(
						newOptionSetProbability[piCount]);
				}
			}
			pMiniStepC = pMiniStepC->pNext();
		}

	// If we move A to directly before pLast
		if (miniStepChoice==intLength-1)
		{
			this->pChain()->insertBefore(pMiniStepA, pLast);
		}
		pMiniStepA->reciprocalRate(rrA);
		pMiniStepA->logOptionSetProbability(losprA);
		pMiniStepA->logChoiceProbability(lcprA);
	}

	this->recordOutcome(*pMiniStepA, accept, MOVE, false);

	delete[] pi;
	delete[] updateProb;
	if (this->pModel()->localML())
	{
		delete[] choiceProbabilityA;
		delete[] optionSetProbabilityA;
		delete[] reciprocalRateA;
	}
	delete[] newChoiceProbability;
	delete[] newOptionSetProbability;	
	delete[] newReciprocalRate;
	
	return accept;
}

/**
 * Implements the MH_permute Metropolis-Hastings step.
 */
bool MLSimulation::permute(int c0)
{
	if (this->pChain()->ministepCount() <= 2)
	{
		return false;
	}

	bool accept = false;

	MiniStep * pMiniStepA = this->pChain()->randomMiniStep();

	while (pMiniStepA == this->pChain()->pLast())
	{
		pMiniStepA = this->pChain()->randomMiniStep();
	}

	vector<MiniStep *> interval;
	MiniStep * pMiniStep = pMiniStepA;

	while ((int) interval.size() < c0 && pMiniStep != this->pChain()->pLast())
	{
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
	}

	if (interval.size() <= 1)
	{
		return false;
	}

	MiniStep * pNextMiniStep = pMiniStep;

	permuteVector(interval);

	this->lthisPermutationLength = interval.size();

	this->setStateBefore(pMiniStepA);
	bool valid = true;
	double sumlprob = 0;
	double sumlprob_new = 0;
	double mu_new = this->pChain()->mu();
	double sigma2_new = this->pChain()->sigma2();
	double * newReciprocalRate = new double[interval.size()];
	double * newOptionSetProbability = new double[interval.size()];
	double * newChoiceProbability = new double[interval.size()];

	for (unsigned i = 0; i < interval.size() && valid; i++)
	{
		pMiniStep = interval[i];
		DependentVariable * pVariable =
			this->lvariables[pMiniStep->variableId()];

		if (!pVariable->validMiniStep(pMiniStep))
		{
			valid = false;
		}
		else
		{
			sumlprob += pMiniStep->logChoiceProbability() +
				pMiniStep->logOptionSetProbability();
			double rrOld = pMiniStep->reciprocalRate();

			if (!this->simpleRates())
			{
				mu_new -= rrOld;
				sigma2_new -= rrOld * rrOld;
			}

			this->calculateRates();
			double rr = 1 / this->grandTotalRate();
			double lospr =
				log(pVariable->rate(pMiniStep->ego()) * rr);
			double lcpr = log(pVariable->probability(pMiniStep));

			sumlprob_new += lospr + lcpr;

			if (!this->simpleRates())
			{
				mu_new += rr;
				sigma2_new += rr * rr;
			}

			pMiniStep->makeChange(pVariable);

			newReciprocalRate[i] = rr;
			newOptionSetProbability[i] = lospr;
			newChoiceProbability[i] = lcpr;
		}
	}

	if (valid)
	{
		double kappaFactor = 1;

		if (!this->simpleRates())
		{
			double sigma2 = this->pChain()->sigma2();
			double mu = this->pChain()->mu();

			kappaFactor = sqrt(sigma2 / sigma2_new) *
				exp((1 - mu) * (1 - mu) / (2 * sigma2) -
					(1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
		}

		this->lproposalProbability =
			kappaFactor * exp(sumlprob_new - sumlprob);

		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}

		accept = nextDouble() < this->lproposalProbability;
		if (accept)
		{
			for (unsigned i = 0; i < interval.size(); i++)
			{
				pMiniStep = interval[i];

				this->pChain()->remove(pMiniStep);

				pMiniStep->reciprocalRate(newReciprocalRate[i]);
				pMiniStep->logOptionSetProbability(
					newOptionSetProbability[i]);
				pMiniStep->logChoiceProbability(newChoiceProbability[i]);
			}

			for (unsigned i = 0; i < interval.size(); i++)
			{
				this->pChain()->insertBefore(interval[i], pNextMiniStep);
			}
		}
	}

	this->recordOutcome(*interval[0], accept, PERMUTE, false);

	delete[] newChoiceProbability;
	delete[] newOptionSetProbability;
	delete[] newReciprocalRate;

	return accept;
}


/**
 * Implements the MH_insPermute Metropolis-Hastings step.
 */
bool MLSimulation::insertPermute(int c0)
{
    if (this->pChain()->ministepCount() <= 1)
	{
		return false;
	}

    bool accept = false;

    MiniStep * pMiniStepA = this->pChain()->randomMiniStep();

    while (pMiniStepA == this->pChain()->pLast())
    {
		pMiniStepA = this->pChain()->randomMiniStep();
    }

    this->setStateBefore(pMiniStepA);
    this->calculateRates();

    DependentVariable * pVariable =
	    this->lvariables[pMiniStepA->variableId()];
    //	double pr2 = 1 - pVariable->rate(pMiniStepA->ego()) / this->grandTotalRate();
    double pr2 = 1;

    pVariable = this->chooseVariable();
    int i = this->chooseActor(pVariable);

    DependentVariable * pVariableInsert = pVariable;

    if (pVariable == this->lvariables[pMiniStepA->variableId()]
		&& i == pMiniStepA->ego())
    {
		return false;
    }

    if (!pVariable->pActorSet()->active(i))
    {
		return false;
    }

    if (pVariable->pData()->upOnly(pVariable->period()) ||
		pVariable ->pData()->downOnly(pVariable->period()))
    {
		return false;
    }

    MiniStep * pLeftMiniStep = pVariable->randomMiniStep(i);

    // 	if (pLeftMiniStep->pOption() == pMiniStepA->pOption())
    // 	{
    // // surely need to delete pLeftMiniStep here too? in fact have quit already if even var and actor
    // are the same as ministep a.
    // 		Rprintf("should have memory leack\n");
    // 		return false;
    // 	}
    if (pLeftMiniStep->diagonal() ||
		!pVariable->validMiniStep(pLeftMiniStep))
    {
		delete pLeftMiniStep;
		return false;
    }

    double rr0 = 1 / this->grandTotalRate();
    double lospr0 = log(pVariable->rate(i) * rr0);
    double lcpr0 = pLeftMiniStep->logChoiceProbability();

    pLeftMiniStep->reciprocalRate(rr0);
    pLeftMiniStep->logOptionSetProbability(lospr0);

    bool misdat = pLeftMiniStep->missing(this->pChain()->period());
    if (misdat && !pLeftMiniStep->missingEnd(this->pChain()->period()))
    {
		delete pLeftMiniStep;
		return false;
    }
    MiniStep * pMiniStepB = this->pChain()->pLast();
    double choiceLength = 1;

    if (!misdat)
    {
		MiniStep * pMiniStepD =
		    this->pChain()->nextMiniStepForOption(*pLeftMiniStep->pOption(),
				pMiniStepA);

		if (!pMiniStepD)
		{
			pMiniStepD = this->pChain()->pLast();
		}

		if (pMiniStepD == pMiniStepA)
		{
			delete pLeftMiniStep;
			return false;
		}

		choiceLength =
		    this->pChain()->intervalLength(pMiniStepA, pMiniStepD) - 1;
		pMiniStepB =
		    this->pChain()->randomMiniStep(pMiniStepA->pNext(), pMiniStepD);
    }

    vector<MiniStep *> interval;
    MiniStep * pMiniStep = pMiniStepA;

    bool validInterval = true;
    while ((int) interval.size() < c0 && pMiniStep != pMiniStepB &&
		validInterval)
    {
		if (!pMiniStep->diagonal())
		{
			for (unsigned i = 0; i < interval.size(); i++)
			{
				if (pMiniStep->pOption() == interval[i]->pOption())
				{
					validInterval = false;
					break;
				}
			}
		}
		if (validInterval)
		{
			interval.push_back(pMiniStep);
			pMiniStep = pMiniStep->pNext();
		}
    }

    permuteVector(interval);

    this->lthisPermutationLength = interval.size();

    // Add the ministeps up to the pMiniStepB to the interval,
    // as the calculations with these ministeps are the same as with
    // the permuted ones.

    while (pMiniStep != pMiniStepB)
    {
		//		Rf_PrintValue(getMiniStepDF(*pMiniStep));
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
    }

    // Variables for localML option.
    int egoInsert = 0;
    int alterInsert = 0;
    NetworkVariable * pNetworkVariableInsert = 0;
    bool * N2vector = 0;
    bool N2 = true;

    if (this->pModel()->localML())
    {
        pNetworkVariableInsert =
			dynamic_cast<NetworkVariable *>(pVariableInsert);
        N2vector = new bool[interval.size()];
        egoInsert = pLeftMiniStep->ego();
        if (pLeftMiniStep->networkMiniStep())
        {
            alterInsert = dynamic_cast<NetworkChange *>(pLeftMiniStep)->alter();
        }
    }

    bool valid = true;
    double sumlprob = 0;
    double sumlprob_new = 0;
    double mu_new = this->pChain()->mu();
    double sigma2_new = this->pChain()->sigma2();
    double * newReciprocalRate = new double[interval.size()];
    double * newOptionSetProbability = new double[interval.size()];
    double * newChoiceProbability = new double[interval.size()];

    // We execute the new chain, do the necessary calculations, and
    // simultaneously check if the new chain is valid.

    // Execute the new ministep that is inserted before pMiniStepA

    sumlprob_new += lcpr0 + lospr0;

    if (!this->simpleRates())
    {
		mu_new += rr0;
		sigma2_new += rr0 * rr0;
    }

    pLeftMiniStep->makeChange(pVariableInsert);

    MiniStep * pRightMiniStep = 0;

    // Execute the existing ministeps (with some of them in a permuted order)

    for (unsigned i = 0; i < interval.size() && valid; i++)
    {
		pMiniStep = interval[i];
		pVariable = this->lvariables[pMiniStep->variableId()];

		if (!pVariable->validMiniStep(pMiniStep))
		{
			valid = false;
			N2 = false;//???
		}
		else
		{
			if (this->pModel()->localML())
			{
				if (i < this->lthisPermutationLength)
				{
					N2 = true;
				} else if (i == this->lthisPermutationLength)
				{
					pRightMiniStep = pLeftMiniStep->createReverseMiniStep();
					pRightMiniStep->makeChange(pVariableInsert);
				}

				if (!(i < this->lthisPermutationLength))
				{
					N2 = this->neighbourhoodChange(pLeftMiniStep, pMiniStep,
						pVariableInsert, pNetworkVariableInsert,
						egoInsert, alterInsert);
				}
			}

			if (!N2)
			{
				pMiniStep->makeChange(pVariable);
			} else
			{
				if ((this->pModel()->localML()) && !(i < this->lthisPermutationLength))
				{
					pLeftMiniStep->makeChange(pVariableInsert);
				}
				sumlprob += pMiniStep->logChoiceProbability() +
					pMiniStep->logOptionSetProbability();

				double rrOld = pMiniStep->reciprocalRate();

				if (!this->simpleRates())
				{
					mu_new -= rrOld;
					sigma2_new -= rrOld * rrOld;
				}

				this->calculateRates();
				double rr = 1 / this->grandTotalRate();
				double lospr =
					log(pVariable->rate(pMiniStep->ego()) * rr);
				double lcpr = log(pVariable->probability(pMiniStep));

				sumlprob_new += lospr + lcpr;

				if (!this->simpleRates())
				{
					mu_new += rr;
					sigma2_new += rr * rr;
				}

				pMiniStep->makeChange(pVariable);

				newReciprocalRate[i] = rr;
				newOptionSetProbability[i] = lospr;
				newChoiceProbability[i] = lcpr;

				if ((this->pModel()->localML()) && !(i < this->lthisPermutationLength))
				{
					pRightMiniStep->makeChange(pVariableInsert);
				}
			}
		}
        if (this->pModel()->localML())
        {
            N2vector[i] = N2;
        }
    }

    // Need to re-execute left step if the interval size is bigger
    // than this->lthisPermutationLength

    if ((interval.size() > this->lthisPermutationLength) &&
		(this->pModel()->localML()))
    {
		pLeftMiniStep->makeChange(pVariableInsert);
    }

    // Execute the canceling ministep of the generated CCP.

    pRightMiniStep = 0;

    if (valid && !misdat)
    {
		pRightMiniStep = pLeftMiniStep->createReverseMiniStep();

		this->calculateRates();
		double rr = 1 / this->grandTotalRate();
		double lospr =
		    log(pVariableInsert->rate(pRightMiniStep->ego()) * rr);

		double lcpr = log(pVariableInsert->probability(pRightMiniStep));

		sumlprob_new += lospr + lcpr;

		if (!this->simpleRates())
		{
			mu_new += rr;
			sigma2_new += rr * rr;
		}

		pRightMiniStep->makeChange(pVariableInsert);
		pRightMiniStep->reciprocalRate(rr);
		pRightMiniStep->logOptionSetProbability(lospr);
		pRightMiniStep->logChoiceProbability(lcpr);
    }

    if (valid)
    {
		double pr1 = 0;

		if (!misdat)
		{
			// Calculate new number of consecutive canceling pairs.
			// Only need to consider those of the type being created as
			// permutation is not allowed to alter the number of others.

			// Network: basic is add 1 if no other option of this type
			// otherwise, 2 unless preceding or succeeding ministeps to A or B
			// are the same.

			// 			int newConsecutiveCancelingPairCount =
			// 				this->pChain()->consecutiveCancelingPairCount() + 1;

			// 			if (this->pChain()->nextMiniStepForOption(
			// 					*(pLeftMiniStep->pOption()), this->pChain()->pFirst()))
			// 			{
			// 				newConsecutiveCancelingPairCount += 1;
			// 				if (pLeftMiniStep->networkMiniStep())
			// 				{
			// 					if (pMiniStepA->pPrevious()->pOption() ==
			// 						pLeftMiniStep->pOption())
			// 					{
			// 						newConsecutiveCancelingPairCount--;
			// 					}
			// 					if (pMiniStepB->pOption() ==
			// 						pLeftMiniStep->pOption())
			// 					{
			// 						newConsecutiveCancelingPairCount--;
			// 					}
			// 				}
			// 				else // behavior miniStep
			// 				{

			// 					BehaviorChange * pPreviousMiniStep =
			// 						dynamic_cast<BehaviorChange *>
			// 						(this->pChain()->nextMiniStepForOption(
			// 							*(pLeftMiniStep->pOption()), pMiniStepA));
			// 					if (pPreviousMiniStep)
			// 					{
			// 						pPreviousMiniStep = dynamic_cast<BehaviorChange *>
			// 							(pPreviousMiniStep->pPreviousWithSameOption());
			// 					}
			// 					BehaviorChange * pNextMiniStep =
			// 						dynamic_cast<BehaviorChange *>
			// 						(this->pChain()->nextMiniStepForOption(
			// 							*(pLeftMiniStep->pOption()), pMiniStepB));
			// 					BehaviorChange * pThisMiniStep =
			// 						dynamic_cast <BehaviorChange *>
			// 						(pLeftMiniStep);
			// 					int d0 = pThisMiniStep->difference();
			// 					int dMinus = 0;
			// 					int dPlus = 0;
			// 					if (pPreviousMiniStep)
			// 					{
			// 						dPlus = pPreviousMiniStep->difference();
			// 					}
			// 					if (pNextMiniStep)
			// 					{
			// 						dMinus = pNextMiniStep->difference();
			// 					}
			// 					if (pPreviousMiniStep == pMiniStepA->pPrevious() ||
			// 						dMinus == d0)
			// 					{
			// 						newConsecutiveCancelingPairCount--;
			// 					}
			// 					if (pNextMiniStep == pMiniStepB ||
			// 						dPlus == d0)
			// 					{
			// 						newConsecutiveCancelingPairCount--;
			// 					}
			// 				}
			// 			}
			// 			Rf_PrintValue(getChainDF(*(this->pChain())));
			// 			this->pChain()->printConsecutiveCancelingPairs();
			this->pChain()->insertBefore(pLeftMiniStep, pMiniStepA);
			this->pChain()->insertBefore(pRightMiniStep, pMiniStepB);
			// 			Rf_PrintValue(getChainDF(*(this->pChain())));
			// 			this->pChain()->printConsecutiveCancelingPairs();
			// 			if (newConsecutiveCancelingPairCount !=
			// 				this->pChain()->consecutiveCancelingPairCount())
			// 			{
			// 				Rprintf("ins diff %d %d \n",newConsecutiveCancelingPairCount,
			// 					this->pChain()->consecutiveCancelingPairCount() );
			// 			Rf_PrintValue(getMiniStepDF(*pMiniStepA));
			//  				Rf_PrintValue(getMiniStepDF(*pMiniStepB));
			//  				Rf_PrintValue(getMiniStepDF(*pLeftMiniStep));
			//  				Rf_PrintValue(getMiniStepDF(*pRightMiniStep));
			//  				Rf_error("time to stop");
			//  			}

			pr1 =
			    (1 -
					this->lmissingNetworkProbability -
					this->lmissingBehaviorProbability) /
			    this->pChain()->consecutiveCancelingPairCount();
			this->pChain()->remove(pLeftMiniStep);
			this->pChain()->remove(pRightMiniStep);
		}
		else
		{
			if (pVariableInsert->networkVariable())
			{
				pr1 = this->lmissingNetworkProbability /
				    (this->pChain()->missingNetworkMiniStepCount() + 1);
			}
			else
			{
				// pVariable is a behavior variable

				pr1 = this->lmissingBehaviorProbability /
				    (this->pChain()->missingBehaviorMiniStepCount() + 1);
			}
		}

		double kappaFactor = 0;

		if (this->simpleRates())
		{
			if (!misdat)
			{
				kappaFactor = 1 /
				    (rr0 * rr0 *
						this->pChain()->ministepCount() *
						(this->pChain()->ministepCount() + 1));
			}
			else
			{
				kappaFactor = 1 / (rr0 * this->pChain()->ministepCount());
			}
		}
		else
		{
			double sigma2 = this->pChain()->sigma2();
			double mu = this->pChain()->mu();

			kappaFactor = sqrt(sigma2 / sigma2_new) *
			    exp((1 - mu) * (1 - mu) / (2 * sigma2) -
				    (1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
		}

		this->lproposalProbability =
		    kappaFactor *
			exp(sumlprob_new - sumlprob) *
		    this->pModel()->deletePermuteProbability() *
			pr1 * pr2 *
			(this->pChain()->ministepCount() - 1) *
			choiceLength /
		    (this->pModel()->insertPermuteProbability() * exp(lospr0 + lcpr0));

		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}

		// 	if(misdat)
		// 	{
		// 		Rf_PrintValue(getMiniStepDF(*pLeftMiniStep));
		// 	}
		accept =  nextDouble() < this->lproposalProbability;

		this->recordOutcome(*pLeftMiniStep,  accept, INSPERM, misdat);
		if (misdat)
		{
			//	Rprintf("mis %d %f \n", accept, this->lproposalProbability);
			//Rf_PrintValue(getMiniStepDF(*pLeftMiniStep));
		}

		if (accept)
		{
			// Change the chain permanently

			for (unsigned i = 0; i < interval.size(); i++)
			{
				pMiniStep = interval[i];
				N2 = true;
				if (this->pModel()->localML())
				{
					N2 = N2vector[i];
				}
				if (N2)
				{
					pMiniStep->reciprocalRate(newReciprocalRate[i]);
					pMiniStep->logOptionSetProbability(
						newOptionSetProbability[i]);
					pMiniStep->logChoiceProbability(newChoiceProbability[i]);
				}
				this->pChain()->remove(pMiniStep);
			}

			this->pChain()->insertBefore(pLeftMiniStep, pMiniStepB);

			for (unsigned i = 0; i < interval.size(); i++)
			{
				this->pChain()->insertBefore(interval[i], pMiniStepB);
			}

			if (!misdat)
			{
				this->pChain()->insertBefore(pRightMiniStep, pMiniStepB);
			}
		}
		else
		{
			// Operation not accepted, so we don't need the generated
			// pair of ministeps.
			delete pLeftMiniStep;
			delete pRightMiniStep;
		}
    }
    else
    {
		delete pLeftMiniStep;
    }
    delete[] newChoiceProbability;
    delete[] newOptionSetProbability;
    delete[] newReciprocalRate;
    if (this->pModel()->localML())
    {
        delete[] N2vector;
    }
    return accept;
}

/**
 * Implements the MH_DelPermute Metropolis-Hastings step
 */
bool MLSimulation::deletePermute(int c0)
{
	bool missdat = FALSE;
    missdat = (nextDouble() <
		this->lmissingNetworkProbability +
		this->lmissingBehaviorProbability);
    MiniStep * pMiniStepA = 0;
    MiniStep * pMiniStepB = 0;
    if (!missdat)
    {
		if (this->pChain()->consecutiveCancelingPairCount() == 0)
		{
			return false;
		}
		pMiniStepA = this->pChain()->randomConsecutiveCancelingPair();
		pMiniStepB = pMiniStepA->pNextWithSameOption();
    }
    else
    {
		if (nextDouble() <
			this->lmissingNetworkProbability /
			(this->lmissingNetworkProbability +
				this->lmissingBehaviorProbability))
		{
			if (this->pChain()->missingNetworkMiniStepCount() == 0)
			{
				return false;
			}

			pMiniStepA = this->pChain()->randomMissingNetworkMiniStep();
			if (pMiniStepA->pNext() == this->pChain()->pLast())
			{
				return false;
			}
		}
		else
		{
			if (this->pChain()->missingBehaviorMiniStepCount() == 0)
			{
				return false;
			}

			pMiniStepA = this->pChain()->randomMissingBehaviorMiniStep();
			if (pMiniStepA->pNext() == this->pChain()->pLast())
			{
				return false;
			}
		}

		pMiniStepB = this->pChain()->pLast();

    }

    // Create the permuted interval of ministeps
    vector<MiniStep *> interval;
    MiniStep * pMiniStep = pMiniStepA->pNext();

    bool validInterval = true;
    while ((int) interval.size() < c0 && pMiniStep != pMiniStepB &&
		validInterval)
    {
		if (!pMiniStep->diagonal())
		{
			for (unsigned i = 0; i < interval.size(); i++)
			{
				if (pMiniStep->pOption() == interval[i]->pOption())
				{
					validInterval = false;
					break;
				}
			}
		}
		if (validInterval)
		{
			interval.push_back(pMiniStep);
			pMiniStep = pMiniStep->pNext();
		}
    }

    permuteVector(interval);

    this->lthisPermutationLength = interval.size();

    // Add the ministeps up to pMiniStepB to the interval,
    // as the calculations with these ministeps are the same as with
    // the permuted ones.

    while (pMiniStep != pMiniStepB)
    {
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
    }

    this->setStateBefore(pMiniStepA);

    bool valid = true;

    double sumlprob = pMiniStepA->logChoiceProbability() +
	    pMiniStepA->logOptionSetProbability();
    double sumlprob_new = 0;
    double rr = pMiniStepA->reciprocalRate();
    double mu_new = this->pChain()->mu() - rr;
    double sigma2_new = this->pChain()->sigma2() - rr * rr;
    double * newReciprocalRate = new double[interval.size()];
    double * newOptionSetProbability = new double[interval.size()];
    double * newChoiceProbability = new double[interval.size()];

    if (!missdat)
    {
		sumlprob += pMiniStepB->logChoiceProbability() +
		    pMiniStepB->logOptionSetProbability();

		if (!this->simpleRates())
		{
			rr = pMiniStepB->reciprocalRate();
			mu_new -= rr;
			sigma2_new -= rr * rr;
		}
    }

    DependentVariable * pVariable =
	    this->lvariables[pMiniStepA->pNext()->variableId()];

    this->calculateRates();
    double pr2 = 1 - pVariable->rate(pMiniStepA->pNext()->ego()) /
	    this->grandTotalRate();
    pr2 = 1;

    // Variables for localML option.

    int egoDelete = 0;
    int alterDelete = 0;
    bool * N2vector = 0;
    bool N2 = true;
    DependentVariable * pVariableDelete = 0;
    NetworkVariable * pNetworkVariableDelete = 0;

    if (this->pModel()->localML())
    {
        egoDelete = pMiniStepA->ego();
        if(pMiniStepA->networkMiniStep())
        {
            alterDelete = dynamic_cast<NetworkChange *>(pMiniStepA)->alter();
        }
        N2vector = new bool[interval.size()];
        pVariableDelete = this->lvariables[pMiniStepA->variableId()];
        pNetworkVariableDelete =
			dynamic_cast<NetworkVariable *>(pVariableDelete);
    }

    // We execute the new chain, do the necessary calculations, and
    // simultaneously check if the new chain is valid.
    for (unsigned i = 0; i < interval.size() && valid; i++)
    {

		pMiniStep = interval[i];
		pVariable = this->lvariables[pMiniStep->variableId()];

		if (!pVariable->validMiniStep(pMiniStep))
		{
			valid = false;
			N2 = false;
		}
		else
		{
			if (this->pModel()->localML())
			{
				if (i < this->lthisPermutationLength)
				{
					N2 = true;
				} else if (i == this->lthisPermutationLength)
				{
					pMiniStepA->makeChange(pVariableDelete);
				}
				if (!(i < this->lthisPermutationLength))
				{
					N2 = this->neighbourhoodChange(pMiniStepA, pMiniStep,
						pVariableDelete,
						pNetworkVariableDelete,
						egoDelete, alterDelete);
				}
			}
			if (!N2)
			{
				pMiniStep->makeChange(pVariable);

			} else
			{
				if ((this->pModel()->localML()) &&
                    !(i < this->lthisPermutationLength))
				{
					pMiniStepB->makeChange(pVariableDelete);

				}
				sumlprob += pMiniStep->logChoiceProbability() +
					pMiniStep->logOptionSetProbability();

				double rrOld = pMiniStep->reciprocalRate();

				if (!this->simpleRates())
				{
					mu_new -= rrOld;
					sigma2_new -= rrOld * rrOld;
				}

				this->calculateRates();
				double rr = 1 / this->grandTotalRate();
				double lospr =
					log(pVariable->rate(pMiniStep->ego()) * rr);
				double lcpr = log(pVariable->probability(pMiniStep));

				sumlprob_new += lospr + lcpr;

				if (!this->simpleRates())
				{
					mu_new += rr;
					sigma2_new += rr * rr;
				}

				pMiniStep->makeChange(pVariable);

				newReciprocalRate[i] = rr;
				newOptionSetProbability[i] = lospr;
				newChoiceProbability[i] = lcpr;

				if ((this->pModel()->localML()) &&
                    !(i < this->lthisPermutationLength))
				{
					pMiniStepA->makeChange(pVariableDelete);
				}
			}

		}
        if (this->pModel()->localML())
        {
            N2vector[i] = N2;
        }
    }


    if ((this->pModel()->localML())&&
		(interval.size() > this->lthisPermutationLength))
    {
        pMiniStepB->makeChange(pVariableDelete);
    }

    bool accept = false;
    if (valid)
    {
		double pr1 = 0;

		if (!missdat)
		{
			// Calculate new number of consecutive canceling pairs.
			// Only need to consider those of the type being deleted as
			// permutation is not allowed to alter the number of others.

			// Network: basic is remove 1 if no other option of this type
			// otherwise, 2 unless preceding or succeeding ministeps to A or B
			// are the same.

			// 			Rf_PrintValue(getChainDF(*(this->pChain())));
			// 			this->pChain()->printConsecutiveCancelingPairs();

			// 			int newConsecutiveCancelingPairCount =
			//  				this->pChain()->consecutiveCancelingPairCount() - 1;

			// 			if (this->pChain()->nextMiniStepForOption(
			// 					*(pMiniStepA->pOption()), this->pChain()->pFirst()))
			// 			{
			// 				newConsecutiveCancelingPairCount -= 1;
			// 				if (pMiniStepA->networkMiniStep())
			// 				{
			// 					if (pMiniStepA->pPrevious()->pOption() ==
			// 						pMiniStepA->pOption())
			// 					{
			// 						newConsecutiveCancelingPairCount++;
			// 					}
			// 					if (pMiniStepB->pNext()->pOption() ==
			// 						pMiniStepB->pOption())
			// 					{
			// 						newConsecutiveCancelingPairCount++;
			// 					}
			// 				}
			// 				else // behavior miniStep
			// 				{

			// 					BehaviorChange * pPreviousMiniStep =
			// 						dynamic_cast<BehaviorChange *>
			// 						(this->pChain()->nextMiniStepForOption(
			// 							*(pMiniStepA->pOption()), pMiniStepA));
			// 					if (pPreviousMiniStep)
			// 					{
			// 						pPreviousMiniStep =
			// 						dynamic_cast<BehaviorChange *>
			// 							(pPreviousMiniStep->pPreviousWithSameOption());
			// 					}
			// 					BehaviorChange * pNextMiniStep =
			// 						dynamic_cast<BehaviorChange *>
			// 						(this->pChain()->nextMiniStepForOption(
			// 							*(pMiniStepB->pOption()), pMiniStepB));
			// 					BehaviorChange * pThisMiniStep =
			// 						dynamic_cast <BehaviorChange *>
			// 						(pMiniStepA);
			// 					int d0 = pThisMiniStep->difference();
			// 					int dMinus = 0;
			// 					int dPlus = 0;
			// 					if (pPreviousMiniStep)
			// 					{
			// 						dPlus = pPreviousMiniStep->difference();
			// 					}
			// 					if (pNextMiniStep)
			// 					{
			// 						dMinus = pNextMiniStep->difference();
			// 					}
			// 					if (pPreviousMiniStep == pMiniStepA->pPrevious() ||
			// 						dMinus == d0)
			// 					{
			// 						newConsecutiveCancelingPairCount++;
			// 					}
			// 					if (pNextMiniStep == pMiniStepB->pNext() ||
			// 						dPlus == d0)
			// 					{
			// 						newConsecutiveCancelingPairCount++;
			// 					}
			// 				}
			// 			}
			// insert and delete
			MiniStep * pAfterA = pMiniStepA->pNext();
			this->pChain()->remove(pMiniStepA);
			MiniStep * pAfterB = pMiniStepB->pNext();
			this->pChain()->remove(pMiniStepB);
			// 			if (newConsecutiveCancelingPairCount !=
			// 				this->pChain()->consecutiveCancelingPairCount())
			// 			{
			// 			this->pChain()->printConsecutiveCancelingPairs();
			// 				Rprintf("diff %d %d \n",newConsecutiveCancelingPairCount,
			// 					this->pChain()->consecutiveCancelingPairCount() );
			// 				Rf_PrintValue(getMiniStepDF(*pMiniStepA));
			// 				Rf_PrintValue(getMiniStepDF(*pMiniStepB));
			// 				Rf_PrintValue(getMiniStepDF(*pAfterA));
			// 				Rf_PrintValue(getMiniStepDF(*pAfterB));
			// 				Rf_error("time to stop");
			// 			}
			this->pChain()->insertBefore(pMiniStepA, pAfterA);
			this->pChain()->insertBefore(pMiniStepB, pAfterB);

			pr1 =
			    (1 -
					this->lmissingNetworkProbability -
					this->lmissingBehaviorProbability) /
			    this->pChain()->consecutiveCancelingPairCount();
		}
		else
		{
			if (pMiniStepA->networkMiniStep())
			{
				pr1 = this->lmissingNetworkProbability /
				    this->pChain()->missingNetworkMiniStepCount();
			}
			else
			{
				// pVariable is a behavior variable

				pr1 = this->lmissingBehaviorProbability /
				    this->pChain()->missingBehaviorMiniStepCount();
			}
		}

		// Calculate choiceLength
		double choiceLength = 1;
		double positionChoice = this->pChain()->ministepCount() - 2;
		if (!missdat)
		{
			MiniStep * pMiniStepD = pMiniStepB->pNextWithSameOption();

			if (!pMiniStepD)
			{
				pMiniStepD = this->pChain()->pLast();
			}

			choiceLength =
			    this->pChain()->intervalLength(pMiniStepA, pMiniStepD) - 3;
			positionChoice --;
		}

		double lpr0 = pMiniStepA->logChoiceProbability() +
		    pMiniStepA->logOptionSetProbability();

		double kappaFactor = 0;

		if (this->simpleRates())
		{
			if (!missdat)
			{
				kappaFactor = rr * rr *
				    (this->pChain()->ministepCount() - 1) *
				    (this->pChain()->ministepCount() - 2);
			}
			else
			{
				kappaFactor = rr * (this->pChain()->ministepCount() - 1);
			}
		}
		else
		{
			double sigma2 = this->pChain()->sigma2();
			double mu = this->pChain()->mu();

			kappaFactor = sqrt(sigma2 / sigma2_new) *
			    exp((1 - mu) * (1 - mu) / (2 * sigma2) -
				    (1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
		}

		this->lproposalProbability =
		    kappaFactor *
			exp(sumlprob_new - sumlprob) *
		    this->pModel()->insertPermuteProbability() *
			exp(lpr0) /
		    (this->pModel()->deletePermuteProbability() *
				pr1 * pr2 *
				// 				(this->pChain()->ministepCount() + 1) *
				positionChoice *
				choiceLength);

		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}

		accept = nextDouble() < this->lproposalProbability;

		this->recordOutcome(*pMiniStepA,  accept, DELPERM,
			missdat);

		if (missdat)
		{
			//Rprintf("dat %d %f \n", accept, this->lproposalProbability);
			//Rf_PrintValue(getMiniStepDF(*pMiniStepA));
		}
		if (accept)
		{
			// Change the chain permanently
			//if (missdat)
			//	Rf_PrintValue(getMiniStepDF(*pMiniStepA));
			this->pChain()->remove(pMiniStepA);
			delete pMiniStepA;

			for (unsigned i = 0; i < interval.size(); i++)
			{
				pMiniStep = interval[i];

				this->pChain()->remove(pMiniStep);

				if (this->pModel()->localML())
				{
					N2 = N2vector[i];
				} else {
					N2 = true;
				}
				if (N2)
				{
					pMiniStep->reciprocalRate(newReciprocalRate[i]);
					pMiniStep->logOptionSetProbability(
						newOptionSetProbability[i]);
					pMiniStep->logChoiceProbability(newChoiceProbability[i]);
				}
			}

			for (unsigned i = 0; i < interval.size(); i++)
			{
				this->pChain()->insertBefore(interval[i], pMiniStepB);
			}

			if (!missdat)
			{
				this->pChain()->remove(pMiniStepB);
				delete pMiniStepB;
			}
		}
    }

    if (this->pModel()->localML())
    {
        delete[] N2vector;
    }
    delete[] newChoiceProbability;
    delete[] newOptionSetProbability;
    delete[] newReciprocalRate;

    return accept;
}


/**
 * Implements the MH_InsMis Metropolis-Hastings step.
 */
bool MLSimulation::insertMissing()
{
	// The numbers in the comments are referencing the specification in
	// Siena_algorithms.pdf.

	// Part A

	// 1.

	if (this->linitialMissingOptions.size() == 0)
	{
		return false;
	}

	// 2.

	const Option * pOption =
		this->linitialMissingOptions[
			nextInt(this->linitialMissingOptions.size())];

	DependentVariable * pVariable =
		this->lvariables[pOption->variableIndex()];
	BehaviorLongitudinalData * pBehaviorData =
		dynamic_cast<BehaviorLongitudinalData *>(pVariable->pData());
	NetworkVariable * pNetworkVariable =
		dynamic_cast<NetworkVariable *>(pVariable);
	BehaviorVariable * pBehaviorVariable =
		dynamic_cast<BehaviorVariable *>(pVariable);

	double pr1 = 1;
	double d0 = 0;

	if (pVariable->behaviorVariable())
	{
		pr1 = 0.5;
		d0 = nextInt(2) * 2 - 1;
	}

	// 3.

	bool reversed = false;

	if (pVariable->behaviorVariable())
	{
		int initialValue =
			this->pChain()->pInitialState()->behaviorValues(pVariable->name())[
				pOption->ego()];
		//double newValue = initialValue + d0;
		double newValue = initialValue - d0;

		if (newValue < pBehaviorData->min() || newValue > pBehaviorData->max())
		{
			pr1 = 1;
			d0 = -d0;
			reversed = true;
		}
	}

	// 4.

	MiniStep * pMiniStepB =
		this->pChain()->firstMiniStepForOption(*pOption);

	if (!pMiniStepB)
	{
		pMiniStepB = this->pChain()->pLast();
	}

	// 5.

	double choiceLength =
		this->pChain()->intervalLength(this->pChain()->pFirst(), pMiniStepB) -
			1;

	// 6.

	MiniStep * pMiniStepA =
		this->pChain()->randomMiniStep(this->pChain()->pFirst()->pNext(),
			pMiniStepB);

	// 7.

	if (pVariable->constrained())
	{
		if (!this->validInsertMissingStep(pOption, d0, pMiniStepA))
		{
			if (pVariable->networkVariable() || reversed)
			{
				return false;
			}
			else
			{
				d0 = -d0;
				pr1 = 1;

				if (!this->validInsertMissingStep(pOption, d0, pMiniStepA))
				{
					return false;
				}
			}
		}
	}

	// Part B

	// 8.

	double sumlprob = 0;
	double sumlprob_new = 0;
	double mu_new = this->pChain()->mu();
	double sigma2_new = this->pChain()->sigma2();

	for (MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
		pMiniStep != pMiniStepA;
		pMiniStep = pMiniStep->pNext())
	{
		sumlprob += pMiniStep->logChoiceProbability() +
			pMiniStep->logOptionSetProbability();
//		Rprintf("insert %f %f %f \n",sumlprob,
//pMiniStep->logChoiceProbability(),
//			pMiniStep->logOptionSetProbability());
		mu_new -= pMiniStep->reciprocalRate();
		sigma2_new -= pMiniStep->reciprocalRate() *
			pMiniStep->reciprocalRate();
	}

	// 9.

	this->setStateBefore(this->pChain()->pFirst()->pNext());
	int initialValue;
	int changedInitialValue;

	if (pVariable->networkVariable())
	{
		initialValue =
			pNetworkVariable->pNetwork()->tieValue(pOption->ego(),
				pOption->alter());
		changedInitialValue = 1 - initialValue;
	}
	else
	{
		initialValue =
			pBehaviorVariable->value(pOption->ego());
		changedInitialValue = initialValue - d0;
	}

	double pr2 =
		pVariable->pData()->observedDistribution(initialValue,
			this->period());
	double pr3 =
		pVariable->pData()->observedDistribution(changedInitialValue,
			this->period());
//	pr2=1;
//	pr3=1;
	// The ministep to be inserted before pMiniStepA
	MiniStep * pNewMiniStep = this->createMiniStep(pOption, d0, false);

	// The dummy ministeps that would change the initial state
	MiniStep * pDummyMiniStep = pNewMiniStep->createReverseMiniStep();

	pDummyMiniStep->makeChange(pVariable);

	// 10.

	int size =
		this->pChain()->intervalLength(this->pChain()->pFirst()->pNext(),
			pMiniStepA) - 1;
	double * newReciprocalRate = new double[size];
	double * newOptionSetProbability = new double[size];
	double * newChoiceProbability = new double[size];

	int i = 0;

	for (MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
		pMiniStep != pMiniStepA;
		pMiniStep = pMiniStep->pNext())
	{
		pVariable = this->lvariables[pMiniStep->variableId()];

		this->calculateRates();
		double rr = 1 / this->grandTotalRate();
		double lospr =
			log(pVariable->rate(pMiniStep->ego()) * rr);
		double lcpr = log(pVariable->probability(pMiniStep));

		sumlprob_new += lospr + lcpr;

//		Rprintf("insert new %f %f %f \n",sumlprob_new,
//			lospr,
//			lcpr);
		if (!this->simpleRates())
		{
			mu_new += rr;
			sigma2_new += rr * rr;
		}

		pMiniStep->makeChange(pVariable);

		newReciprocalRate[i] = rr;
		newOptionSetProbability[i] = lospr;
		newChoiceProbability[i] = lcpr;

		i++;
	}

	// 11.

	pVariable = this->lvariables[pNewMiniStep->variableId()];
	this->calculateRates();
	double rr0 = 1 / this->grandTotalRate();
	double lospr0 =
		log(pVariable->rate(pNewMiniStep->ego()) * rr0);
	double lcpr0 = log(pVariable->probability(pNewMiniStep));

	sumlprob_new += lospr0 + lcpr0;

	if (!this->simpleRates())
	{
		mu_new += rr0;
		sigma2_new += rr0 * rr0;
	}

	pNewMiniStep->reciprocalRate(rr0);
	pNewMiniStep->logChoiceProbability(lcpr0);
	pNewMiniStep->logOptionSetProbability(lospr0);

	// 12.

	double kappaFactor;
	double mu = this->pChain()->mu();
	double sigma2 = this->pChain()->sigma2();

	if (this->simpleRates())
	{
		kappaFactor = 1 / (rr0 * this->pChain()->ministepCount());
	}
	else
	{
		kappaFactor =
			sqrt(sigma2 / sigma2_new) *
				exp((1 - mu) * (1 - mu) / (2 * sigma2) -
					(1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
	}

	// 13.
//	Rprintf("insert %f %f %f %f %f %f %f %f %f \n", kappaFactor,
	//	sumlprob_new, sumlprob,
	//	this->pModel()->insertRandomMissingProbability(),
	//	pr1, pr3, pr2, choiceLength,
	//	this->pModel()->deleteRandomMissingProbability());

	this->lproposalProbability =
		kappaFactor *
			exp(sumlprob_new - sumlprob) *
			this->pModel()->deleteRandomMissingProbability() *
			choiceLength *
			pr3 /
		(this->pModel()->insertRandomMissingProbability() * pr1 * pr2);

	if (this->lproposalProbability > 1)
	{
		this->lproposalProbability = 1;
	}

	// Part C

	// 14.

	bool accept = nextDouble() < this->lproposalProbability;

	this->recordOutcome(*pNewMiniStep, accept, INSMISS, false);

	// 15.
	//Rprintf("ins %d %f\n", accept, this->lproposalProbability);
	//Rf_PrintValue(getMiniStepDF(*pNewMiniStep));

	if (accept)
	{
		// 	Rf_PrintValue(getMiniStepDF(*pNewMiniStep));
		this->pChain()->changeInitialState(pDummyMiniStep);

		int i = 0;

		for (MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
			pMiniStep != pMiniStepA;
			pMiniStep = pMiniStep->pNext())
		{
			pMiniStep->logChoiceProbability(newChoiceProbability[i]);
			pMiniStep->logOptionSetProbability(newOptionSetProbability[i]);
			pMiniStep->reciprocalRate(newReciprocalRate[i]);
			i++;
		}

		this->pChain()->insertBefore(pNewMiniStep, pMiniStepA);
	}
	else
	{
		// The ministep was not accepted, so we don't need it.
		delete pNewMiniStep;
	}

	delete[] newReciprocalRate;
	delete[] newOptionSetProbability;
	delete[] newChoiceProbability;
	delete pDummyMiniStep;

	return accept;
}


/**
 * Implements the MH_DelMis Metropolis-Hastings step.
 */
bool MLSimulation::deleteMissing()
{
	// The numbers in the comments are referencing the specification in
	// Siena_algorithms.pdf.

	// Part A

	// 1.

	if (this->linitialMissingOptions.size() == 0)
	{
		return false;
	}

	// 2.

	const Option * pOption =
		this->linitialMissingOptions[
			nextInt(this->linitialMissingOptions.size())];

	// 3.

	MiniStep * pMiniStepA =
		this->pChain()->firstMiniStepForOption(*pOption);

	if (!pMiniStepA)
	{
		return false;
	}

	int d0 = 0;

	if (pMiniStepA->behaviorMiniStep())
	{
		BehaviorChange * pBehaviorChange =
			dynamic_cast<BehaviorChange *>(pMiniStepA);
		d0 = pBehaviorChange->difference();
	}

	// 4.

	MiniStep * pMiniStepB = pMiniStepA->pNextWithSameOption();

	if (!pMiniStepB)
	{
		pMiniStepB = this->pChain()->pLast();
	}

	double choiceLength =
		this->pChain()->intervalLength(this->pChain()->pFirst(), pMiniStepB) - 2;

	// 5.

	DependentVariable * pVariable =
		this->lvariables[pOption->variableIndex()];
	BehaviorLongitudinalData * pBehaviorData =
		dynamic_cast<BehaviorLongitudinalData *>(pVariable->pData());

	double pr1 = 0.5;

	if (pVariable->behaviorVariable())
	{
		int initialValue =
			this->pChain()->pInitialState()->behaviorValues(pVariable->name())[
				pOption->ego()];
		double testValue = initialValue + 2 * d0;

		if (testValue < pBehaviorData->min() ||
			testValue > pBehaviorData->max())
		{
			pr1 = 1;
		}
	}
	else
	{
		pr1 = 1.0;
	}
	// 6.

	if (pVariable->constrained())
	{
		if (!this->validDeleteMissingStep(pMiniStepA, false))
		{
			return false;
		}

		if (pr1 == 0.5 && !this->validDeleteMissingStep(pMiniStepA, true))
		{
			pr1 = 1;
		}
	}

	// Part B

	// 7.

	double sumlprob = 0;
	double sumlprob_new = 0;
	double mu_new = this->pChain()->mu();
	double sigma2_new = this->pChain()->sigma2();

	for (MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
		pMiniStep != pMiniStepA->pNext();
		pMiniStep = pMiniStep->pNext())
	{
		sumlprob += pMiniStep->logChoiceProbability() +
			pMiniStep->logOptionSetProbability();
//		Rprintf("delete %f %f %f \n",sumlprob,
//pMiniStep->logChoiceProbability(),
//			pMiniStep->logOptionSetProbability());
		mu_new -= pMiniStep->reciprocalRate();
		sigma2_new -= pMiniStep->reciprocalRate() *
			pMiniStep->reciprocalRate();
	}

	// 8.

	NetworkVariable * pNetworkVariable =
		dynamic_cast<NetworkVariable *>(pVariable);
	BehaviorVariable * pBehaviorVariable =
		dynamic_cast<BehaviorVariable *>(pVariable);

	this->resetVariables();
	int initialValue;
	int changedInitialValue;

	if (pVariable->networkVariable())
	{
		initialValue =
			pNetworkVariable->pNetwork()->tieValue(pOption->ego(),
				pOption->alter());
		changedInitialValue = 1 - initialValue;
	}
	else
	{
		initialValue =
			pBehaviorVariable->value(pOption->ego());
		changedInitialValue = initialValue + d0;
	}

	double pr2 =
		pVariable->pData()->observedDistribution(initialValue,
			this->period());
	double pr3 =
		pVariable->pData()->observedDistribution(changedInitialValue,
			this->period());
//	pr2=1;
//	pr3=1;
	this->calculateRates();
	double rr0 = 1 / this->grandTotalRate();
//	double lospr0 =
//		log(pVariable->rate(pOption->ego()) * rr0);
//	double lcpr0 = log(pVariable->probability(pMiniStepA));

	//sumlprob_new += lospr0 + lcpr0;

	if (!this->simpleRates())
	{
		mu_new += rr0;
		sigma2_new += rr0 * rr0;
	}

	pMiniStepA->makeChange(pVariable);

	// 9.

	int size =
		this->pChain()->intervalLength(this->pChain()->pFirst()->pNext(),
			pMiniStepA) - 1;
	double * newReciprocalRate = new double[size];
	double * newOptionSetProbability = new double[size];
	double * newChoiceProbability = new double[size];

	int i = 0;

	for (MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
		pMiniStep != pMiniStepA;
		pMiniStep = pMiniStep->pNext())
	{
		pVariable = this->lvariables[pMiniStep->variableId()];

		this->calculateRates();
		double rr = 1 / this->grandTotalRate();
		double lospr =
			log(pVariable->rate(pMiniStep->ego()) * rr);
		double lcpr = log(pVariable->probability(pMiniStep));

		sumlprob_new += lospr + lcpr;
//	Rprintf("delete new %f %f %f \n",sumlprob_new,
//lospr, lcpr);

		if (!this->simpleRates())
		{
			mu_new += rr;
			sigma2_new += rr * rr;
		}

		pMiniStep->makeChange(pVariable);

		newReciprocalRate[i] = rr;
		newOptionSetProbability[i] = lospr;
		newChoiceProbability[i] = lcpr;

		i++;
	}

	// 10.

	double kappaFactor;

	if (this->simpleRates())
	{
		kappaFactor = rr0 * (this->pChain()->ministepCount() - 1);
	}
	else
	{
		double mu = this->pChain()->mu();
		double sigma2 = this->pChain()->sigma2();

		kappaFactor =
			sqrt(sigma2 / sigma2_new) *
				exp((1 - mu) * (1 - mu) / (2 * sigma2) -
					(1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
	}

	// 11.
//	Rprintf("delete %f %f %f %f %f %f %f %f %f \n", kappaFactor,
//		sumlprob_new, sumlprob,
//		this->pModel()->insertRandomMissingProbability(),
//		pr1, pr3, pr2, choiceLength,
//		this->pModel()->deleteRandomMissingProbability());
	this->lproposalProbability =
		kappaFactor *
			exp(sumlprob_new - sumlprob) *
			this->pModel()->insertRandomMissingProbability() *
			pr1 *
			pr3 /
		(this->pModel()->deleteRandomMissingProbability() *
			choiceLength *
			pr2);

	if (this->lproposalProbability > 1)
	{
		this->lproposalProbability = 1;
	}

	// Part C

	// 12.

	bool accept = nextDouble() < this->lproposalProbability;

	// 13.
	this->recordOutcome(*pMiniStepA, accept, DELMISS, false);
	//Rprintf("del %d %f\n", accept, this->lproposalProbability);
	//Rf_PrintValue(getMiniStepDF(*pMiniStepA));

	if (accept)
	{
		//	Rf_PrintValue(getMiniStepDF(*pMiniStepA));
		this->pChain()->changeInitialState(pMiniStepA);

		int i = 0;

		for (MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
			pMiniStep != pMiniStepA;
			pMiniStep = pMiniStep->pNext())
		{
			pMiniStep->logChoiceProbability(newChoiceProbability[i]);
			pMiniStep->logOptionSetProbability(newOptionSetProbability[i]);
			pMiniStep->reciprocalRate(newReciprocalRate[i]);
			i++;
		}

		this->pChain()->remove(pMiniStepA);
		delete pMiniStepA;
	}

	delete[] newReciprocalRate;
	delete[] newOptionSetProbability;
	delete[] newChoiceProbability;

	return accept;
}


/**
 * Returns if inserting the (w0, i0, j0, r0, d0) right before the given
 * mini step pMiniStepA and changing the initial state by (w0, i0, j0, r0, -d0)
 * is valid with respect to various constraints. The option (w0, i0, j0, r0)
 * is given by the parameter pOption.
 */
bool MLSimulation::validInsertMissingStep(const Option  * pOption,
	int d0,
	const MiniStep * pMiniStepA)
{
	bool rc = true;

	// Initialize the variables
	this->resetVariables();

	DependentVariable * pVariable = this->lvariables[pOption->variableIndex()];

	// The ministep to be inserted before pMiniStepA
	MiniStep * pNewMiniStep = this->createMiniStep(pOption, d0, false);

	// The dummy ministeps that would change the initial state
	MiniStep * pDummyMiniStep = pNewMiniStep->createReverseMiniStep();

	// Check if changing the initial state according to the dummy ministep
	// is valid. Note that the up-only and down-only conditions should not
	// be tested, as this is not an actual change, but rather just a
	// convenient way of changing the initial state.

	if (pVariable->validMiniStep(pDummyMiniStep, false))
	{
		pDummyMiniStep->makeChange(pVariable);
	}
	else
	{
		rc = false;
	}

	// Test the validity of existing ministeps up to pMiniStepA and
	// execute them.

	for (MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
		pMiniStep != pMiniStepA && rc;
		pMiniStep = pMiniStep->pNext())
	{
		pVariable = this->lvariables[pMiniStep->variableId()];

		if (pVariable->validMiniStep(pMiniStep))
		{
			pMiniStep->makeChange(pVariable);
		}
		else
		{
			rc = false;
		}
	}

	// Finally, test the new ministep to be inserted before pMiniStepA

	if (rc)
	{
		pVariable = this->lvariables[pNewMiniStep->variableId()];
		rc = pVariable->validMiniStep(pNewMiniStep);
	}

	// Free the memory

	delete pDummyMiniStep;
	delete pNewMiniStep;

	return rc;
}


/**
 * Returns if deleting the given ministep and changing the initial state by
 * the same ministep to compensate the removal is valid with respect to various
 * constraints. If the parameter <code>applyTwice</code> is <code>true</code>,
 * then the given ministep is applied two times initially and instead of being
 * deleted, it is replaced with the reverse ministep. So this apply-twice mode
 * doubles the effect of a single delete-missing step.
 */
bool MLSimulation::validDeleteMissingStep(MiniStep * pMiniStepA,
	bool applyTwice)
{
	bool rc = true;

	// Initialize the variables
	this->resetVariables();

	DependentVariable * pVariable =
		this->lvariables[pMiniStepA->variableId()];

	// Check if changing the initial state according to the given ministep
	// is valid. Note that the up-only and down-only conditions should not
	// be tested, as this is not an actual change, but rather just a
	// convenient way of changing the initial state.

	if (pVariable->validMiniStep(pMiniStepA, false))
	{
		pMiniStepA->makeChange(pVariable);
	}
	else
	{
		rc = false;
	}

	if (applyTwice)
	{
		if (pVariable->validMiniStep(pMiniStepA, false))
		{
			pMiniStepA->makeChange(pVariable);
		}
		else
		{
			rc = false;
		}
	}

	// Test the validity of existing ministeps up to the given ministep and
	// execute them.

	for (MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
		pMiniStep != pMiniStepA && rc;
		pMiniStep = pMiniStep->pNext())
	{
		pVariable = this->lvariables[pMiniStep->variableId()];

		if (pVariable->validMiniStep(pMiniStep))
		{
			pMiniStep->makeChange(pVariable);
		}
		else
		{
			rc = false;
		}
	}

	if (rc && applyTwice)
	{
		MiniStep * pReverseMiniStep = pMiniStepA->createReverseMiniStep();
		pVariable = this->lvariables[pReverseMiniStep->variableId()];
		rc = pVariable->validMiniStep(pReverseMiniStep);
		delete pReverseMiniStep;
	}

	return rc;
}


/**
 * Creates and returns a new mini step for the given option and difference
 * parameter (the latter for behavior changes only).
 */
MiniStep * MLSimulation::createMiniStep(const Option * pOption,
	int difference, bool diagonal) const
{
	MiniStep * pMiniStep = 0;
	DependentVariable * pVariable = this->lvariables[pOption->variableIndex()];

	if (pVariable->networkVariable())
	{
		pMiniStep =
			new NetworkChange(
				dynamic_cast<NetworkLongitudinalData *>(pVariable->pData()),
				pOption->ego(),
				pOption->alter(), diagonal);
	}
	else
	{
		pMiniStep =
			new BehaviorChange(
				dynamic_cast<BehaviorLongitudinalData *>(pVariable->pData()),
				pOption->ego(),
				difference);
	}

	return pMiniStep;
}


/**
 * increments the number of aborted steps
 * for the given steptype for this variable
 */
void MLSimulation::incrementAborts(int stepType)
{
	this->laborts[stepType]++;
}

// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the proposal probability of the last Metropolis-Hastings step.
 */
double MLSimulation::proposalProbability() const
{
	return this->lproposalProbability;
}



/**
 * Stores the missing network probability (prmin in siena_Algorithms).
 */
void MLSimulation::missingNetworkProbability(double probability)
{
	this->lmissingNetworkProbability = probability;
}


/**
 * Returns the missing network probability (prmin in siena_Algorithms).
 */
double MLSimulation::missingNetworkProbability() const
{
	return this->lmissingNetworkProbability;
}


/**
 * Stores the missing behavior probability (prmib in siena_Algorithms).
 */
void MLSimulation::missingBehaviorProbability(double probability)
{
	this->lmissingBehaviorProbability = probability;
}


/**
 * Returns the missing behavior probability (prmib in siena_Algorithms).
 */
double MLSimulation::missingBehaviorProbability() const
{
	return this->lmissingBehaviorProbability;
}


/**
 * Returns the number of acceptances for this steptype.
 */
//int MLSimulation::acceptances(int stepType) const
//{
//	return this->lacceptances[stepType];
//}

/**
 * Returns the number of rejections for this steptype.
 */
//int MLSimulation::rejections(int stepType) const
//{
//	return this->lrejections[stepType];
//}
/**
 * Returns the number of aborted steps for this steptype.
 */
int MLSimulation::aborts(int stepType) const
{
	return this->laborts[stepType];
}

/**
 * Updates the permutation length for this period.
 */

void MLSimulation::updateCurrentPermutationLength(bool accept)
{
	unsigned permutationLength = this->lcurrentPermutationLength;
	if (this->lthisPermutationLength == permutationLength)
	{
		double minvalue = this->pModel()->minimumPermutationLength();
		double maxvalue = this->pModel()->maximumPermutationLength();
		if (accept)
		{
			this->lcurrentPermutationLength += 0.5;
			if (this->lcurrentPermutationLength > maxvalue)
			{
				this->lcurrentPermutationLength = maxvalue;
			}
		}
		else
		{
			this->lcurrentPermutationLength -= 0.5;
			if (this->lcurrentPermutationLength < minvalue)
			{
				this->lcurrentPermutationLength = minvalue;
			}
		}
	}
}
/**
 * Returns the permutation length for this period.
 */

double MLSimulation::currentPermutationLength() const
{
	return  this->lcurrentPermutationLength;
}
/**
 * Stores the permutation length for this period.
 */

void MLSimulation::currentPermutationLength(double value)
{
	this->lcurrentPermutationLength = value;
}
/**
 * Generates a vector of ministeps representing the differences
 * between the current state at the end of the period and the observation.
 */
void MLSimulation::createEndStateDifferences()
{
	this->pChain()->clearEndStateDifferences();
	const Data * pData = this->pData();
	int period = this->period();

//	Rf_PrintValue(getChainDF(*this->pChain()));
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
		const DependentVariable * pDependentVariable =
			this->pVariable(pVariableData->name());

		if (pNetworkData)
		{
			const Network * pNetwork1 = pNetworkData->pNetwork(period + 1);
			const NetworkVariable * pNetworkVariable =
				dynamic_cast<const NetworkVariable *>(pDependentVariable);

			const Network * pNetwork2 =	pNetworkVariable->pNetwork();

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
								period + 1)
							//	|| !pNetworkData->structural(i, iter1.actor(),
							//	period + 1)
							)
						{
							NetworkChange * pMiniStep =
								new NetworkChange(pNetworkData,
									i,
									iter1.actor(), false);
							//	if (!pNetworkData->missing(i, iter1.actor(),
							//		period + 1))
							//{
							//	Rf_PrintValue(getMiniStepDF(*pMiniStep));
								//}
							this->pChain()->addEndStateDifference(pMiniStep);

							iter1.next();
							//	Rf_PrintValue(getMiniStepDF(
							//		*this->lendStateDifferences.back()));
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
								period + 1)
							//	|| !pNetworkData->structural(i, iter2.actor(),
							// period + 1)
							)
						{
							this->pChain()->addEndStateDifference(
							//	this->lendStateDifferences.push_back(
								new NetworkChange(pNetworkData,
									i,
									iter2.actor(), false));
							//	if (!pNetworkData->missing(i, iter2.actor(),
							//			period + 1))
							//{
							//	Rf_PrintValue(getMiniStepDF(
							//			*this->lendStateDifferences.back()));
								//}
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
			const BehaviorVariable * pBehaviorVariable =
				dynamic_cast<const BehaviorVariable *>(pDependentVariable);

			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				int delta = pBehaviorVariable->value(i)
					- pBehaviorData->value(period + 1, i);
				int singleChange = 1;

				if (delta < 0)
				{
					delta = -delta;
					singleChange = -1;
				}

				for (int j = 0; j < delta; j++)
				{
					if (!pBehaviorData->structural(period + 1, j)
						//|| !pBehaviorData->structural(period, j + 1)
						)

					{
						this->pChain()->addEndStateDifference(
							//this->lendStateDifferences.push_back(
							new BehaviorChange(pBehaviorData,
								i,
								singleChange));
						//			Rprintf(" %d %d in beh\n", i, singleChange);
						//			Rf_PrintValue(getMiniStepDF(
						//					*this->lendStateDifferences.back()));
					}
					else
					{
						// create step in structural subchain?
					}
				}
			}
		}
	}
	// Rprintf("xx %d %d %d end create end state diff\n",
	// 	this->lendStateDifferences.size(), period, this->pChain()->ministepCount() - 1);
}


/**
 * Calculates whether toggling a ministep affects the small neighbourhood of the
 * ego of another ministep.
*/
bool MLSimulation::smallNeighbourhoodChange(MiniStep * pMiniStep1,
												MiniStep * pMiniStep2,
												DependentVariable * pVariable,
												NetworkVariable * pNetworkVariable,
												int ego1, int alter1)
{
	//// Won't work for behaviour and multiplex networks.
	bool addToInterval = false;
		int ego = pMiniStep2->ego();
	int alter = dynamic_cast<NetworkChange *>(pMiniStep2)->alter();
		if ((ego==ego1) || (ego==alter1) || (alter==ego1) || (alter==alter1))
	{
		addToInterval = true;
	}
	return addToInterval;
}


/**
 * Calculates whether toggling a ministep affects the neighbourhood of the ego
 * of another ministep.
 */
bool MLSimulation::neighbourhoodChange(MiniStep * pMiniStep1,
	MiniStep * pMiniStep2,
	DependentVariable * pVariable,
	NetworkVariable * pNetworkVariable,
	int ego1, int alter1)
{
	bool addToInterval = false;

	if (this->lvariables[pMiniStep2->variableId()]->pActorSet() ==
		pVariable->pActorSet())
	{
		NetworkVariable * pNetworkVariable2 = 0;
		int myties = 0;
		bool net = false;
		int ego = pMiniStep2->ego();

		if (this->lvariables[pMiniStep2->variableId()]== pVariable)
		{
			// Loop over all networks with same actorset
			if (this->lvariables.size() > 1)
			{
				for (unsigned k = 0;k < this->lvariables.size(); k++)
				{
					if (this->lvariables[k]->pActorSet()==
						pVariable->pActorSet())
					{
						pNetworkVariable2 = dynamic_cast<NetworkVariable *>(
							this->lvariables[k]);

						if (pNetworkVariable2)
						{
							if (pNetworkVariable2->pNetwork()->
								tieValue(ego,ego1))
							{
								net = true;
							}
						}
					}
				}
			}
		}
		if ((pMiniStep2->networkMiniStep())&&(pMiniStep1->networkMiniStep()))
		{
			if (pNetworkVariable->pNetwork()->tieValue(ego1,ego) ||
				pNetworkVariable->pNetwork()->tieValue(alter1,ego) ||
				pNetworkVariable->pNetwork()->tieValue(ego,ego1) ||
				pNetworkVariable->pNetwork()->tieValue(ego,alter1))
			{myties = 1;} else {myties = 0;}

			if (this->lvariables[pMiniStep2->variableId()]== pVariable)
			{
				// Both ministeps change the same network

				if ( (ego==ego1) || (ego==alter1) ||
					(myties==1) || net)
				{
					addToInterval = true;
				}
			} else {
				// Multiplex: ministeps change different network with same
				// actor sets

				if ( (ego==ego1) || (ego==alter1) ||
					(myties==1))
				{
					addToInterval = true;
				}
			}

		} else if (pMiniStep1->networkMiniStep()) {
			// Ministep1 is a network change; ministep2 a behaviour change

			if ((ego==ego1) || (ego==alter1) ||
				(pNetworkVariable->pNetwork()->tieValue(ego,alter1)))
			{
				addToInterval = true;
			}
		} else if (pMiniStep2->networkMiniStep()) {
			// Ministep1 is a behaviour change, ministep2 a network change
			addToInterval = true;
		} else if (this->lvariables[pMiniStep2->variableId()]== pVariable){
			// Both ministeps changes the same behaviour
			if ((ego==ego1) || net)
			{
				addToInterval = true;
			}
		} else {
			// Ministeps change different behaviours
			addToInterval = true;
		}
	}

	return addToInterval;
}

/**
 * Returns the simulated data at the end of the period.
 */
void MLSimulation::gotoLastState()
{
	this->setStateBefore(this->pChain()->pLast());
}

}
