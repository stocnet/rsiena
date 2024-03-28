/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkVariable.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkVariable class.
 *****************************************************************************/
#include <algorithm>
#include <vector>
#include <cmath>
#include <R_ext/Print.h>
#include <R_ext/Arith.h>
#include <Rinternals.h>
#include "NetworkVariable.h"
#include "network/NetworkUtils.h"
#include "utils/Utils.h"
#include "utils/Random.h"
#include "data/ActorSet.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "network/iterators/UnionTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/SimulationActorSet.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "model/effects/NetworkEffect.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/ml/NetworkChange.h"
#include "model/filters/PermittedChangeFilter.h"
#include "model/ml/Chain.h"
#include "model/settings/Setting.h"

#include <algorithm>
#include <vector>
#include <cmath>

#include <R_ext/Error.h>
#include <R_ext/Print.h>
#include <R_ext/Arith.h>
#include <Rinternals.h>

using namespace std;

namespace siena
{

SEXP getMiniStepDF(const MiniStep& miniStep);

// ----------------------------------------------------------------------------
// Section: Construction
// ----------------------------------------------------------------------------

/**
 * Creates a new network variable for the given observed data.
 * @param pSimulation the owner of this variable
 */
NetworkVariable::NetworkVariable(NetworkLongitudinalData * pData,
	EpochSimulation * pSimulation) :
		DependentVariable(pData->name(),
			pData->pActorSet(),
			pSimulation)
{
	this->lpData = pData;
	this->lpSenders = pSimulation->pSimulationActorSet(pData->pSenders());
	this->lpReceivers = pSimulation->pSimulationActorSet(pData->pReceivers());
	this->lpNetwork = 0;
	this->lactiveStructuralTieCount = new int[this->n()];
	int numberOfAlters;
	this->loneMode = pData->oneModeNetwork();
	if (this->loneMode)
	{
		this->lpNetwork = new OneModeNetwork(this->n(), false);
		numberOfAlters = this->m();
		this->lpermitted = new bool[this->m()];
		this->levaluationEffectContribution = new double * [numberOfAlters];
		this->lendowmentEffectContribution = new double * [numberOfAlters];
		this->lcreationEffectContribution = new double * [numberOfAlters];
		this->lprobabilities = new double[numberOfAlters];
	}
	else
	{
		this->lpNetwork = new Network(this->n(), this->m());
		numberOfAlters = this->m() + 1;
		// then numberOfAlters really is a misnomer
		this->lpermitted = new bool[numberOfAlters];
		this->levaluationEffectContribution = new double * [numberOfAlters];
		this->lendowmentEffectContribution = new double * [numberOfAlters];
		this->lcreationEffectContribution = new double * [numberOfAlters];
		this->lprobabilities = new double[numberOfAlters];
	}
	for (int i = 0; i < numberSettings(); i++) {
		lsettings[i]->initSetting(lpNetwork);
	}
	this->lsymmetricEvaluationEffectContribution = new double * [2];
	this->lsymmetricEndowmentEffectContribution = new double * [2];
	this->lsymmetricCreationEffectContribution = new double * [2];

	for (int i = 0; i < numberOfAlters; i++)
	{
		this->levaluationEffectContribution[i] =
			new double[pSimulation->pModel()->
				rEvaluationEffects(pData->name()).size()];
		this->lendowmentEffectContribution[i] =
			new double[pSimulation->pModel()->
				rEndowmentEffects(pData->name()).size()];
		this->lcreationEffectContribution[i] =
			new double[pSimulation->pModel()->
				rCreationEffects(pData->name()).size()];
	}
	for (int i = 0; i < 2; i++)
	{
		this->lsymmetricEvaluationEffectContribution[i] =
			new double[pSimulation->pModel()->
				rEvaluationEffects(pData->name()).size()];
		this->lsymmetricEndowmentEffectContribution[i] =
			new double[pSimulation->pModel()->
				rEndowmentEffects(pData->name()).size()];
		this->lsymmetricCreationEffectContribution[i] =
			new double[pSimulation->pModel()->
				rCreationEffects(pData->name()).size()];
	}

	this->lpNetworkCache =
		pSimulation->pCache()->pNetworkCache(this->lpNetwork);

	this->lalter = 0;

	this->lnetworkModelType = NetworkModelType(pData->modelType());
	this->lnetworkModelTypeDoubleStep = false;
	this->lnetworkDoubleStepProb = 0;

	if (this->loneMode)
	{
		if (this->lnetworkModelType == DOUBLESTEP25)
		{
			this->lnetworkModelTypeDoubleStep = true;
			this->lnetworkDoubleStepProb = 0.25;
		}
		if (this->lnetworkModelType == DOUBLESTEP50)
		{
			this->lnetworkModelTypeDoubleStep = true;
			this->lnetworkDoubleStepProb = 0.50;
		}
		if (this->lnetworkModelType == DOUBLESTEP75)
		{
			this->lnetworkModelTypeDoubleStep = true;
			this->lnetworkDoubleStepProb = 0.75;
		}
		if (this->lnetworkModelType == DOUBLESTEP100)
		{
			this->lnetworkModelTypeDoubleStep = true;
			this->lnetworkDoubleStepProb = 1.00;
		}
	}
}


/**
 * Deallocates this variable object.
 */
NetworkVariable::~NetworkVariable()
{
	for (int i = 0; i < numberSettings(); i++) {
		lsettings[i]->terminateSetting(lpNetwork);
	}
	delete this->lpNetwork;

	delete[] this->lactiveStructuralTieCount;
	delete[] this->lpermitted;
	delete[] this->lprobabilities;

	// Delete arrays of contributions
	int numberOfAlters;
	if (this->oneModeNetwork())
	{
		numberOfAlters = this->m();
	}
	else
	{
		numberOfAlters = this->m() + 1;
	}
	for (int i = 0; i < numberOfAlters; i++)
	{
		delete[] this->levaluationEffectContribution[i];
		delete[] this->lendowmentEffectContribution[i];
		delete[] this->lcreationEffectContribution[i];
	}

	for (int i = 0; i < 2; i++)
	{
		delete[] this->lsymmetricEvaluationEffectContribution[i];
		delete[] this->lsymmetricEndowmentEffectContribution[i];
		delete[] this->lsymmetricCreationEffectContribution[i];
	}

	delete[] this->levaluationEffectContribution;
	delete[] this->lendowmentEffectContribution;
	delete[] this->lcreationEffectContribution;

	delete[] this->lsymmetricEvaluationEffectContribution;
	delete[] this->lsymmetricEndowmentEffectContribution;
	delete[] this->lsymmetricCreationEffectContribution;

	this->lsymmetricEvaluationEffectContribution = 0;
	this->lsymmetricEndowmentEffectContribution = 0;
	this->lsymmetricCreationEffectContribution = 0;
	this->levaluationEffectContribution = 0;
	this->lendowmentEffectContribution = 0;
	this->lcreationEffectContribution = 0;
	this->lpData = 0;
	this->lpNetwork = 0;
	this->lactiveStructuralTieCount = 0;
	this->lpermitted = 0;
	this->lprobabilities = 0;

	// no need to delete lpChangeContribution since this is
	// handled by the MiniStep

	deallocateVector(this->lpermittedChangeFilters);
}


/**
 * Adds the given filter of permissible changes to the filters of this
 * variable. This variable becomes the owner of the filter, which means that
 * the filter will be deleted as soon as this variable is deleted.
 */
void NetworkVariable::addPermittedChangeFilter(PermittedChangeFilter * pFilter)
{
	this->lpermittedChangeFilters.push_back(pFilter);
}



// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the set of actors acting as tie senders.
 */
const SimulationActorSet * NetworkVariable::pSenders() const
{
	return this->lpSenders;
}


/**
 * Returns the set of actors acting as tie receivers.
 */
const SimulationActorSet * NetworkVariable::pReceivers() const
{
	return this->lpReceivers;
}


/**
 * Returns the second dimension of this variable, namely, how many values
 * correspond to each actor. This number equals the number of receivers for
 * network variables.
 */
int NetworkVariable::m() const
{
	return this->lpReceivers->n();
}


/**
 * Indicates if this is a one-mode network, an attribute of the network
 */
bool NetworkVariable::oneModeNetwork() const
{
	return this->loneMode;
}


/**
 * Returns the longitudinal data object this variable is based on.
 */
LongitudinalData * NetworkVariable::pData() const
{
	return this->lpData;
}

/**
 * Sets the model type.
 */
void NetworkVariable::networkModelType(int type)
{
	this->lnetworkModelType = NetworkModelType(type);
}

/**
 * Returns the model type.
 */
NetworkModelType NetworkVariable::networkModelType() const
{
	return this->lnetworkModelType;
}

/**
 * Returns whether the model type is one of the symmetric type b models.
 */

bool NetworkVariable::networkModelTypeB() const
{
	return this->lnetworkModelType == BFORCE ||
		this->lnetworkModelType == BAGREE || this->lnetworkModelType == BJOINT;
}

/**
 * Returns whether the model type is one of the double step models.
 */

bool NetworkVariable::networkModelTypeDoubleStep() const
{
	return this->lnetworkModelTypeDoubleStep;
}


/**
 * Returns the probability for the DOUBLESTEP model.
 */
double NetworkVariable::networkDoubleStepProb() const
{
	return this->lnetworkDoubleStepProb;
}


// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Initializes this variable as of the beginning of the given period.
 */
void NetworkVariable::initialize(int period)
{
	DependentVariable::initialize(period);

	// Copy the respective observation

	if (this->oneModeNetwork())
	{
		OneModeNetwork * pNetwork =
			(OneModeNetwork *) this->lpNetwork;
		const OneModeNetwork * pObservedNetwork =
			(const OneModeNetwork *) this->lpData->pNetwork(period);

		// Use the copy assignment operator
		(*pNetwork) = (*pObservedNetwork);
	}
	else
	{
		// Use the copy assignment operator
		(*this->lpNetwork) = (*this->lpData->pNetwork(period));
	}

	// Initialize the counters with all structural ties, including those to
	// inactive actors.

	for (int i = 0; i < this->n(); i++)
	{
		this->lactiveStructuralTieCount[i] =
			this->lpData->structuralTieCount(i, period);
	}

	// Now subtract the structural ties to initially inactive actors.

	for (int i = 0; i < this->m(); i++)
	{
		if (!this->pReceivers()->active(i))
		{
			for (IncidentTieIterator iter =
					this->lpData->pStructuralTieNetwork(period)->inTies(i, "nwva");
				iter.valid();
				iter.next())
			{
				this->lactiveStructuralTieCount[iter.actor()]--;
			}
		}
	}

}


/**
 * Returns the current state of the network.
 */
Network * NetworkVariable::pNetwork() const
{
	return this->lpNetwork;
}


// ----------------------------------------------------------------------------
// Section: Composition change
// ----------------------------------------------------------------------------

/**
 * Updates the current network and other variables when an actor becomes
 * active.
 */
void NetworkVariable::actOnJoiner(const SimulationActorSet * pActorSet,
	int actor)
{
	DependentVariable::actOnJoiner(pActorSet, actor);

	const Network * pStartNetwork = this->lpData->pNetwork(this->period());

	if (pActorSet == this->pSenders())
	{
		// Activate the ties to active receivers according to the
		// initial observation of the period.

		for (IncidentTieIterator iter = pStartNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			if (this->pReceivers()->active(iter.actor()))
			{
				this->lpNetwork->setTieValue(actor,
					iter.actor(),
					iter.value());
			}
		}

		// The rates need to be recalculated.
		this->invalidateRates();
	}

	if (pActorSet == this->pReceivers())
	{
		// Activate the ties from active senders according to the
		// initial observation of the period.

		for (IncidentTieIterator iter = pStartNetwork->inTies(actor, "nwvb");
			iter.valid();
			iter.next())
		{
			if (this->pSenders()->active(iter.actor()))
			{
				this->lpNetwork->setTieValue(iter.actor(),
					actor,
					iter.value());
			}
		}

		// Update the numbers of structural ties to active actors, as one actor
		// becomes active now.

		for (IncidentTieIterator iter =
				this->lpData->pStructuralTieNetwork(this->period())->inTies(
					actor, "nwvc");
			iter.valid();
			iter.next())
		{
			this->lactiveStructuralTieCount[iter.actor()]++;
		}

		// The rates need to be recalculated.
		this->invalidateRates();
	}
}


/**
 * Updates the current network and other variables when an actor becomes
 * inactive.
 */
void NetworkVariable::actOnLeaver(const SimulationActorSet * pActorSet,
	int actor)
{
	DependentVariable::actOnLeaver(pActorSet, actor);

	if (pActorSet == this->pSenders())
	{
		// Remove the ties from the given actor.
		this->lpNetwork->clearOutTies(actor);

		// The rates need to be recalculated.
		this->invalidateRates();
	}

	if (pActorSet == this->pReceivers())
	{
		// Remove the ties to the given actor.
		this->lpNetwork->clearInTies(actor);

		// Update the numbers of structural ties to active actors, as one actor
		// becomes inactive now.

		for (IncidentTieIterator iter =
				this->lpData->pStructuralTieNetwork(this->period())->inTies(
					actor, "nwvd");
			iter.valid();
			iter.next())
		{
			this->lactiveStructuralTieCount[iter.actor()]--;
		}

		// The rates need to be recalculated.
		this->invalidateRates();
	}
}


/**
 * Sets leavers values back to the value at the start of the simulation.
 *
 */
void NetworkVariable::setLeaverBack(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->pSenders())
	{
		// Reset ties from the given actor to values at start

		for (int i = 0; i < this->m(); i++)
		{
			if (i != actor)
			{
				this->lpNetwork->setTieValue(actor,
					i,
					this->lpData->tieValue(actor, i, this->period()));
			}
		}
	}

	if (pActorSet == this->pReceivers())
	{
		for (int i = 0; i < this->n(); i++)
		{
			if (i != actor)
			{
				this->lpNetwork->setTieValue(i,
					actor,
					this->lpData->tieValue(i, actor, this->period()));
			}
		}
	}
}


// ----------------------------------------------------------------------------
// Section: Changing the network
// ----------------------------------------------------------------------------

/**
 * Returns if the given actor can change the current state of this variable,
 * namely, it is active and at least one of the tie variables to active
 * receivers is not structurally determined.
 */
bool NetworkVariable::canMakeChange(int actor) const
{
	bool rc = DependentVariable::canMakeChange(actor);

	if (rc)
	{
		int activeAlterCount =
			this->lpReceivers->activeActorCount();

		if (this->oneModeNetwork())
		{
			// No loops are possible in one mode networks
			activeAlterCount--;
		}

		rc &= this->lpSenders->active(actor) &&
			this->lactiveStructuralTieCount[actor] < activeAlterCount;
	}

	return rc;
}


/**
 * Simulates a change of the network according to the choice of the given
 * actor. First, the actor chooses the alter based on the evaluation and
 * endowment functions, and then the tie to the selected alter is flipped.
 * Used for MoM.
 */
void NetworkVariable::makeChange(int actor)
{
	int m;
	this->lego = actor;
	bool accept = true;
	int alter;

	this->successfulChange(true);

	if (this->symmetric() && this->networkModelTypeB())
	{
		if (this->calculateModelTypeBProbabilities())
		{
			accept = nextDouble() < this->lsymmetricProbability;

			if (this->pSimulation()->pModel()->needScores())
			{
				this->accumulateSymmetricModelScores(this->lalter, accept);
			}
		}
		else
		{
			this->successfulChange(false);
			return;
		}
		alter = this->lalter;
// this->lalter is determined in calculateModelTypeBProbabilities()
	}
	else
	{
		this->calculateTieFlipProbabilities();

		if (this->oneModeNetwork())
		{
			m = this->m();
		}
		else
		{
			m = this->m() + 1;
		}
		if (this->stepType() != -1)
		{
			if (!lsettings[stepType()]->validate(lpNetwork)) {
				successfulChange(false);
				lsettings[stepType()]->terminateSetting(ego());
				return;
			}
		}

		alter = nextIntWithProbabilities(m,	this->lprobabilities);

		if (this->lnetworkModelTypeDoubleStep)
		{
			this->lalter = alter;
		}

//		if (this->stepType() > 1) {
//			if (this->lpNetworkCache->outTieExists(alter)) {
//				this->successfulChange(false);
//				return;
//			}
//		}
		// Siena 3 checks in the diagonal case, so I do too temporarily.
		//if (alter != actor && !this->lpNetworkCache->outTieExists(alter) &&
		//	this->pSimulation()->pModel()->modelType() == AAGREE)
		if (this->symmetric() &&
			this->networkModelType() == AAGREE &&
			!this->lpNetworkCache->outTieExists(alter))
		{
			this->checkAlterAgreement(alter);
			double value = nextDouble();
			accept = value < this->lsymmetricProbability;

			if (this->pSimulation()->pModel()->needScores())
			{
				this->addAlterAgreementScores(accept);
			}
		}

		if (this->pSimulation()->pModel()->needScores())
		{
			this->accumulateScores(alter);
		}
		if (this->pSimulation()->pModel()->needDerivatives())
		{
			this->accumulateDerivatives(); // ABC
		}
	}
//	 NB  the probabilities in the reported chain are probably wrong for !accept
	if (this->pSimulation()->pModel()->needChain())
	{
		// add ministep to chain
		MiniStep * pMiniStep;
		if (accept)
		{
			pMiniStep =
				new NetworkChange(this->lpData, actor, alter,
					this->diagonalMiniStep(actor, alter));
		}
		else
		{
			pMiniStep =
				new NetworkChange(this->lpData, actor, alter, true);
		}
		if (this->pSimulation()->pModel()->needChangeContributions())
		{
			pMiniStep->changeContributions(lpChangeContribution);
		}
		this->pSimulation()->pChain()->insertBefore(pMiniStep,
			this->pSimulation()->pChain()->pLast());
		if (!this->symmetric() || !this->networkModelTypeB())
		{
			pMiniStep->logChoiceProbability(log(this->lprobabilities[alter]));
			if (this->symmetric() &&
				this->networkModelType() == AAGREE)
			{
				pMiniStep->logChoiceProbability(pMiniStep->
					logChoiceProbability() + log(this->lsymmetricProbability));
			}
		}
		else
		{
			double probability = this->lsymmetricProbability;
			if (!accept)
			{
				probability = 1 - probability;
			}
			pMiniStep->logChoiceProbability(log(this->lalterProbability) +
				log(probability));
		}
	}
	// Make a change if we have a real alter (other than the ego or
	// the dummy for bipartite networks)

	if (accept && ((!this->oneModeNetwork() && alter < this->m()) ||
			(this->oneModeNetwork() && this->lego != alter)))
	{
		int currentValue = this->lpNetwork->tieValue(this->lego, alter);

		// Update the distance from the observed data at the beginning of the
		// period. Ties missing at any of the endpoints of the period
		// don't contribute to the distance

		// If network is symmetric, changes are in steps of two, otherwise 1
		int change = 1;
		if (this->oneModeNetwork())
		{
			if ((dynamic_cast<const OneModeNetworkLongitudinalData *>
					(this->lpData))->symmetric())
			{
				change = 2;
			}
		}
		if (!this->lpData->missing(this->lego, alter, this->period()) &&
			!this->lpData->missing(this->lego, alter, this->period() + 1))
		{
			if (this->lpData->tieValue(this->lego, alter, this->period()) ==
				currentValue)
			{
				this->simulatedDistance(this->simulatedDistance() + change);
			}
			else
			{
				this->simulatedDistance(this->simulatedDistance() - change);
			}
		}

		this->lpNetwork->setTieValue(this->lego, alter, 1 - currentValue);

		if (this->oneModeNetwork())
		{
			const OneModeNetworkLongitudinalData * pData =
				dynamic_cast<const OneModeNetworkLongitudinalData *>(
					this->pData());

			if (pData->symmetric())
			{
				this->lpNetwork->setTieValue(alter,
					this->lego,
					1 - currentValue);
			}
		}
	}
	if (stepType() != -1) {
		lsettings[stepType()]->terminateSetting(ego());
	}
}


/**
 * This method does some preprocessing to speed up subsequent queries regarding
 * the specified (usually current) ego.
 */
void NetworkVariable::preprocessEgo(int ego)
{
	// Let the effects do their preprocessing.

	this->preprocessEgo(this->pEvaluationFunction(), ego);
	this->preprocessEgo(this->pEndowmentFunction(), ego);
	this->preprocessEgo(this->pCreationFunction(), ego);
}


/**
 * This method does some preprocessing for each effect of the given
 * function to speed up subsequent queries regarding
 * the specified (usually current) ego.
 */
void NetworkVariable::preprocessEgo(const Function * pFunction, int ego)
{
	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		NetworkEffect * pEffect =
			(NetworkEffect *) pFunction->rEffects()[i];
		pEffect->preprocessEgo(ego);
	}
}


/**
 * Determines the set of actors that are allowed to act as alters in the next
 * tie flip, and stores this information in <code>lpermitted</code>.
 */
void NetworkVariable::calculatePermissibleChanges()
{
	NetworkLongitudinalData * pData = (NetworkLongitudinalData *) this->pData();

	int m = this->m();
	for (int i = 0; i < m; ++i)
	{
		this->lpermitted[i] = false;
	}

	Setting* curSetting = 0;
	ITieIterator* iter = 0;
	if (this->stepType() != -1)
	{
		curSetting = lsettings[stepType()];
		iter = curSetting->getSteps();
		m = curSetting->getSize();
	}

  // Test each alter if a tie flip to that alter is permitted according to
  // upOnly/downOnly flags and the max degree.
	for (int ii = 0; ii < m; ii++)
	{
		int i = ii;
		if (this->stepType() > -1) {
			if (!iter->valid()) {
				Rf_error( "size of iterator != size setting");
			}
			i = iter->actor();
			iter->next();
		}

		if (this->lpNetworkCache->outTieExists(i))
		{
			this->lpermitted[i] = !pData->upOnly(this->period());
		}
		else if (i != this->lego || !this->oneModeNetwork())
		{
			this->lpermitted[i] =
				// for comparability with siena3 comment this out
				//	this->pSimulation()->active(this->pReceivers(), i) &&
				(!pData->downOnly(this->period())) &&
				this->lpNetwork->outDegree(this->lego) < pData->maxDegree();
		}
		else
		{
			// It is okay to not make any change at all
			this->lpermitted[i] = true;
		}
	}
	if (iter != 0) {
		delete iter;
	}

	// Prohibit the change of structural ties

	for (IncidentTieIterator iter =
			 pData->pStructuralTieNetwork(this->period())->outTies(this->lego);
		 iter.valid();
		 iter.next())
	{
		this->lpermitted[iter.actor()] = false;
	}

	// Run the filters that may forbid some more changes

	for (unsigned i = 0; i < this->lpermittedChangeFilters.size(); i++)
	{
		PermittedChangeFilter * pFilter = this->lpermittedChangeFilters[i];
		//	Rprintf(" about to filter %d %d\n ", i, this->lpermitted[10]);
		pFilter->filterPermittedChanges(this->lego, this->lpermitted);
		//Rprintf(" filtered %d %d\n ", i, this->lpermitted[10]);
	}

	// ensure that no change is a permitted change
	if (this->oneModeNetwork())
	{
		this->lpermitted[lego] = true;
	} else {
		this->lpermitted[m] = true;
	}

}


/**
 * For each alter, this method calculates the contribution of each evaluation,
 * endowment, and tie creation effect if a tie from the ego to this alter was
 * flipped.
 * These contributions are stored in arrays
 * <code>levaluationEffectContribution</code>,
 * <code>lendowmentEffectContribution</code>, and
 * <code>lcreationEffectContribution</code>.
 */
void NetworkVariable::calculateTieFlipContributions()
{
	int evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	int endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();
	int creationEffectCount = this->pCreationFunction()->rEffects().size();
	const vector<Effect *> & rEvaluationEffects =
		this->pEvaluationFunction()->rEffects();
	const vector<Effect *> & rEndowmentEffects =
		this->pEndowmentFunction()->rEffects();
	const vector<Effect *> & rCreationEffects =
		this->pCreationFunction()->rEffects();
	bool twoModeNetwork = !this->oneModeNetwork();

	int m = this->m();
	int alter = 0;
	Setting* curSetting = 0;
	ITieIterator* permIter = 0;
	if (this->stepType() != -1)
	{
		curSetting = lsettings[stepType()];
		permIter = curSetting->getPermittedSteps();
		m = curSetting->getPermittedSize();
	}
	for (int alteri = 0; alteri < m; alteri++)
	{
		alter = alteri;

		if (this->stepType() != -1)
		{
			if (!permIter->valid()) {
				Rf_error("permitted iter length != settings permitted size");
			}
			alter = permIter->actor();
			permIter->next();
		}

		// alter = ego for one-mode networks means no change.
		// alter = m for two-mode networks means no change.
		// No change, zero contribution; this is done at the end
		// of calculateTieFlipContributions.

		if (this->lpermitted[alter] && (twoModeNetwork || alter != this->lego))
		// The second condition seems superfluous,
		// as (alter == this->lego) is dealt with later.
		// TODO: check if this can safely be dropped.
		{
			for (int i = 0; i < evaluationEffectCount; i++)
			{
				NetworkEffect * pEffect =
					(NetworkEffect *) rEvaluationEffects[i];
				double contribution = pEffect->calculateContribution(alter);

				// Tie withdrawals contribute in the opposite way

				if (this->lpNetworkCache->outTieExists(alter))
				{
					contribution = -contribution;
				}

				this->levaluationEffectContribution[alter][i] = contribution;
			}
		}
		else
		{
			if (!this->lpermitted[alter])
				{
				for (int i = 0; i < evaluationEffectCount; i++)
				{
				this->levaluationEffectContribution[alter][i] = R_NaN;
				}
			}
		}

		// The endowment effects have non-zero contributions on tie
		// withdrawals only
		// TODO: Isn't there && (twoModeNetwork || alter != this->lego) missing?
		if (this->lpNetworkCache->outTieExists(alter) &&
			this->lpermitted[alter])
		{
			for (int i = 0; i < endowmentEffectCount; i++)
			{
				NetworkEffect * pEffect =
					(NetworkEffect *) rEndowmentEffects[i];
				this->lendowmentEffectContribution[alter][i] =
					-pEffect->calculateContribution(alter);
			}
		}
		else
		{
			for (int i = 0; i < endowmentEffectCount; i++)
			{
				if (!this->lpermitted[alter])
				{
					this->lendowmentEffectContribution[alter][i] = R_NaN;
				}
				else
				{
					this->lendowmentEffectContribution[alter][i] = 0;
				}
			}
		}

		// The tie creation effects have non-zero contributions on tie
		// creation only

		if (!this->lpNetworkCache->outTieExists(alter) &&
			this->lpermitted[alter] &&
			(twoModeNetwork || alter != this->lego))
		{
			for (int i = 0; i < creationEffectCount; i++)
			{
				NetworkEffect * pEffect =
					(NetworkEffect *) rCreationEffects[i];
				double contribution =
					pEffect->calculateContribution(alter);

				this->lcreationEffectContribution[alter][i] = contribution;
			}
		}
		else
		{
			for (int i = 0; i < creationEffectCount; i++)
			{
				if (!this->lpermitted[alter])
				{
					this->lcreationEffectContribution[alter][i] = R_NaN;
				}
				else
				{
					this->lcreationEffectContribution[alter][i] = 0;
				}
			}
		}
	}
	if (permIter != 0) {
		delete permIter;
	}

	// need to initialise the no-effect option of alter = this->m() for
	// twomode networks

	if (twoModeNetwork)
	{
		for (int i = 0; i < evaluationEffectCount; i++)
		{
			this->levaluationEffectContribution[m][i] = 0;
		}

		for (int i = 0; i < endowmentEffectCount; i++)
		{
			this->lendowmentEffectContribution[m][i] = 0;
		}

		for (int i = 0; i < creationEffectCount; i++)
		{
			this->lcreationEffectContribution[m][i] = 0;
		}
	}
	else
	{
		for (int i = 0; i < evaluationEffectCount; i++)
		{
			this->levaluationEffectContribution[this->lego][i] = 0;
		}

		for (int i = 0; i < endowmentEffectCount; i++)
		{
			this->lendowmentEffectContribution[this->lego][i] = 0;
		}

		for (int i = 0; i < creationEffectCount; i++)
		{
			this->lcreationEffectContribution[this->lego][i] = 0;
		}
	}

}


/**
 * If Settings model, chooses the setting;
 * calculates the probability of each actor
 * for being chosen as alter for the next tie flip.
 */
void NetworkVariable::calculateTieFlipProbabilities()
{
	NetworkLongitudinalData * pData = (NetworkLongitudinalData *) this->pData();
	// If a settings model, decide what to do. Pre-store probabilities
	this->getStepType();
//	if (this->stepType() > 1)
//	{
//		// clear settings vector
//		this->lsetting->clear();
//		// find the covariate and make a vector of potential alters
//		NetworkLongitudinalData * pNetworkData =
//			dynamic_cast<NetworkLongitudinalData *>(this->pData());
//		string settingName = pNetworkData->rSettingNames()[this->stepType()];
//		ConstantDyadicCovariate * pConstantCovariate =
//			this->pSimulation()->pData()->pConstantDyadicCovariate(settingName);
//	   ChangingDyadicCovariate * pChangingCovariate =
//			this->pSimulation()->pData()->pChangingDyadicCovariate(settingName);
//		map<int, double> row;
//		if (pConstantCovariate) // TODO: What is this?
//		{
//			row = pConstantCovariate->rRowValues(this->lego);
//		}
//		else if (pChangingCovariate)
//		{
//			row = pChangingCovariate->rRowValues(this->lego, this->period());
//		}
//		else
//		{
//			Rf_error("setting not found");
//		}
//		bool needEgo = true;
//		if (row.size() == 0)
//		{
//			this->lsetting->push_back(this->lego);
//			needEgo = false;
//		}
//
//		for (map<int, double>::iterator iter=row.begin();
//			 iter != row.end(); iter++ )
//		{
//			if (iter->first > this->lego && needEgo)
//			{
//				this->lsetting->push_back(this->lego);
//				needEgo = false;
//			}
//			this->lsetting->push_back(iter->first);
//		}
//		if (needEgo)
//		{
//			this->lsetting->push_back(this->lego);
//		}
//	}
//	// TODO: correct this->lsetting for primary setting
//	if (this->stepType() == 1)
//	{
//	{
//		delete this->lsetting;
//		this->lsetting = primarySetting(this->pNetwork(), this->lego);
//	}
	if (stepType() > -1)
	{
		initializeSetting();
		this->lpNetworkCache->stepTypeSet(stepType());
	}
	this->preprocessEgo(this->lego);
	this->calculatePermissibleChanges();
//	int m = this->m();
//
//	if (!this->oneModeNetwork())
//	{
//		m = this->m() + 1;
//	}
	// TODO: check correctness
	if (this->stepType() != -1)
	{
		lsettings[stepType()]->initPermittedSteps(lpermitted);
		// if there is only one permissable change (ego) return
		if (!lsettings[stepType()]->validate(lpNetwork))
		{
			return;
		}
	}
//		//  check we can find an ego in the setting
//		int numberPermitted = 0;
//		for (int i = 0; i < m; i++)
//		{
//			int ii = i;
//			if (this->stepType() > 0)
//			{
//				ii = (*this->lsetting)[i];
//			}
//			if (this->lpermitted[ii])
//			{
//				numberPermitted++;
//			}
//		}
//		int fastAlter = this->lego;
//
//		if (numberPermitted > 1)
//		{
//			while(fastAlter == this->lego && this->lpermitted[fastAlter])
//			{
//				fastAlter = nextInt(m);  // random number
//				if (this->stepType()  > 1)
//				{
//					fastAlter = (*this->lsetting)[fastAlter];
//				}
//			}
//		}
//		this->lsetting->clear();
//		this->lsetting->push_back(min(this->lego, fastAlter));
//		this->lsetting->push_back(max(this->lego, fastAlter));
//	}

	this->calculateTieFlipContributions();

	int evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	int endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();
	int creationEffectCount = this->pCreationFunction()->rEffects().size();

	if (this->pSimulation()->pModel()->needChangeContributions()) {
		this->lpChangeContribution =
				new map<const EffectInfo *, vector<double> >();
// used to be  >[evaluationEffectCount+endowmentEffectCount+creationEffectCount];
		for (int i = 0; i < evaluationEffectCount; i++) {
			vector<double> vec(this->m(),0);
			this->lpChangeContribution->insert(make_pair(
						this->pEvaluationFunction()->rEffects()[i]->pEffectInfo(), vec));
		}
		for (int i = 0; i < endowmentEffectCount; i++)
		{
			vector<double> vec(this->m(),0);
			this->lpChangeContribution->insert(make_pair(
						this->pEndowmentFunction()->rEffects()[i]->pEffectInfo(), vec));
		}
		for (int i = 0; i < creationEffectCount; i++)
		{
			vector<double> vec(this->m(),0);
			this->lpChangeContribution->insert(make_pair(
						this->pCreationFunction()->rEffects()[i]->pEffectInfo(), vec));
		}
	}

	double total = 0;
	double maxValue = 0; // the maximum never can be less than 0
					// because there always is the no-change option
	int m = this->m(); // not this->m()+1 for two-mode network;
					// this is handled separately below
//	double primaryOffset = 0;

	Setting* curSetting = 0;
	ITieIterator* permIter = 0;
	if (this->stepType() != -1)
	{
		curSetting = lsettings[stepType()];
		permIter = curSetting->getPermittedSteps();
		m = curSetting->getPermittedSize();
		if (this->stepType() == 1)
		{
			int egoOutDegree = 0;
			egoOutDegree = this->lpNetwork->outDegree(this->lego);
			if (egoOutDegree > m)
			{
					Rf_error("outdegree > primary setting size");
			}
//			else if (egoOutDegree < m)
//			{
//				primaryOffset = -std::log(m - egoOutDegree);
//			}
		}
	}

	// calculate all contributions then subtract the largest to remove overflow
	int alter = 0;
	int sumPermitted = 0;
	for (int alteri = 0; alteri < m; alteri++)
	{
		alter = alteri;
		if (curSetting)
		{
			if (!permIter->valid())
			{
				Rf_error( "permIter size differs from setting size");
			}
			alter = permIter->actor();
			permIter->next();
		}
		if (this->lpermitted[alter])
		{
			// Calculate the total contribution of all effects

			double contribution = 0;

			for (int i = 0; i < evaluationEffectCount; i++)
			{
				Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
				contribution += pEffect->parameter()
					* this->levaluationEffectContribution[alter][i];
				if (this->pSimulation()->pModel()->needChangeContributions())
				{
					(*this->lpChangeContribution)[pEffect->pEffectInfo()].at(alter) =
					this->levaluationEffectContribution[alter][i];
				}
			}

			if (this->lpNetworkCache->outTieExists(alter))
			{
				for (int i = 0; i < endowmentEffectCount; i++)
				{
					Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];
					contribution += pEffect->parameter()
						* this->lendowmentEffectContribution[alter][i];
					if (this->pSimulation()->pModel()->needChangeContributions())
					{
						(* this->lpChangeContribution)[pEffect->pEffectInfo()].at(alter) =
							this->lendowmentEffectContribution[alter][i];
					}
				}
			}
			else
			{
				for (int i = 0; i < creationEffectCount; i++)
				{
					Effect * pEffect = this->pCreationFunction()->rEffects()[i];
					contribution += pEffect->parameter()
						* this->lcreationEffectContribution[alter][i];
					if (this->pSimulation()->pModel()->needChangeContributions())
					{
						(* this->lpChangeContribution)[pEffect->pEffectInfo()].at(alter) =
							this->lcreationEffectContribution[alter][i];
					}
				}
//				for settings model:
				if ((this->stepType() == 0) || (this->stepType() >= 2))
				{
					contribution += pData->universalOffset();
				}
//				else if (this->stepType() == 1)
//				{
//					contribution += primaryOffset;
//				}
			}

			// The selection probability is the exponential of the total
			// contribution.
			this->lprobabilities[alter] = contribution;
			//	Rprintf("p1 %d %f \n", alter, contribution);
			sumPermitted++;
		}
		else
		{
			this->lprobabilities[alter] = R_NegInf;
		}
		maxValue = max(maxValue, this->lprobabilities[alter]);
	}

	if (permIter != 0)
	{
		permIter->reset();
	}

	for (int alteri = 0; alteri < m; alteri++) {
		alter = alteri;

		if (this->stepType() != -1)
		{
			if (!permIter->valid())
			{
				Rf_error( "permitted iter length != settings permitted size");
			}
			alter = permIter->actor();
			permIter->next();
		}

		if (this->lpermitted[alter])
		{
			this->lprobabilities[alter] -= 	maxValue;
			this->lprobabilities[alter] = exp(this->lprobabilities[alter]);
			total += this->lprobabilities[alter];
		}
		else
		{
			this->lprobabilities[alter] = 0;
		}
	}

	// reset m in case of settings
	m = this->m();
	if (!this->oneModeNetwork())
	{
		if (sumPermitted >= 2)
		{
			this->lprobabilities[m] = exp(-maxValue); // no-change option
			total += this->lprobabilities[m];
		}
		else
		{
			this->lprobabilities[m] = 1; // only the no-change option exists
			total = 1;
		}
		m++; // Watch out: this is meant only for the next normalization, keep
		     // this in mind if further statements are added to this procedure.
	}

	// delete iter
	if (permIter != 0)
	{
		permIter->reset();
	}

	// Normalize
	if (total > 0)
	{
		// this loop requires that permIter processes its entries in ascending order (which it does)
		for (int alter = 0; alter < m; alter++)
		{
			if (permIter != 0) // settings model
			{
				if (permIter->valid())
				{
					if (alter < permIter->actor())
					{
						lprobabilities[alter] = 0;
					}
					else
					{
						lprobabilities[alter] /= total;
						permIter->next();
					}
				}
				else
				{
					lprobabilities[alter] = 0;
				}
			} // else non-settings model
			else if (lpermitted[alter])
			{
				this->lprobabilities[alter] /= total;
			}
			else
			{
				this->lprobabilities[alter] = 0;
			}
		}
	}
	else
	{
		Rprintf("total = %f\n", total);
		Rprintf("this actor = %d\n", (this->lego + 1) );
		Rprintf("this period = %d\n", (this->period() + 1) );
		// counting starts at 0
		Rf_error("total probability non-positive");
	}

	// delete iter
	if (permIter != 0) {
		delete permIter;
	}
}

/**
 * Updates the scores for effects according
 * to the current step in the simulation.
 */
void NetworkVariable::accumulateScores(int alter) const
{
	int m = this->m();
	int sumPermitted = 0;
	Setting* curSetting = 0;
	ITieIterator* permIter = 0;

	if (this->stepType() != -1)
	{
		curSetting = lsettings[stepType()];
		m = curSetting->getPermittedSize();
		sumPermitted = m;
		permIter = curSetting->getPermittedSteps();
	} else
	{
		if (!this->oneModeNetwork()) {
		m++;
	}
		if (alter >= m) {
			Rprintf("this->n = %d this->m = %d m = %d alter = %d \n", this->n(),
					this->m(), m, alter);
		Rf_error("alter too large");
	}
	for (int h = 0; h < m; h++)
	{
		if (this->lpermitted[h]) sumPermitted++;
	}
	}
	if (sumPermitted <= 0)
	{
		Rf_error("nothing was permitted");
	}
	else if (sumPermitted >= 2)
		// if sumPermitted == 1, no contribution to scores
	{
		for (unsigned i = 0;
			i < this->pEvaluationFunction()->rEffects().size();
			i++)
		{
			Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
			double score = this->levaluationEffectContribution[alter][i];
			if (R_IsNaN(score))
			{
				Rprintf("R_IsNaN error: i = %d ego = %d alter = %d m = %d\n",
					i, this->lego, alter, m);
				Rf_error("nan score 41");
			}
			if (curSetting) {
				permIter->reset();
			}
			// subtract sum over all permitted
			int j = 0;
			for (int alteri = 0; alteri < m; alteri++)
			{
				j = alteri;
				if (permIter)
				{
					if (!permIter->valid())
					{
						Rf_error("iterator not valid");
					}
					j = permIter->actor();
					permIter->next();
				}
				if (this->lpermitted[j])
				{
					score -=
						this->levaluationEffectContribution[j][i] *
						this->lprobabilities[j];
				}
				if (R_IsNaN(score))
				{
					Rprintf("R_IsNaN error: i = %d ego = %d alter = %d j = %d m = %d\n",
							i, this->lego, alter, j, m);
					//						Rprintf("lpermitted: \n");
					//						for (int hh = 0; hh < m; hh++)
					//						{
					//							if (this->lpermitted[hh])
					//							{
					//								Rprintf("%d 1", hh);
					//							}
					//							else
					//							{
					//								Rprintf("%d 0", hh);
					//							}
					//							Rprintf("\n");
					//						}
					Rprintf("R_IsNaN error: this->levaluationEffectContribution[j][i] = %f\n",
							this->levaluationEffectContribution[j][i]);
					Rprintf("R_IsNaN Rf_error: this->lprobabilities[j] = %f\n",
							this->lprobabilities[j]);
					Rf_error("nan score 1");
				}
			}
			if (R_IsNaN(this->pSimulation()->score(pEffect->pEffectInfo())))
			{
				Rprintf("R_IsNaN error: i = %d ego = %d alter = %d m = %d\n",
					i, this->lego, alter, m);
					Rf_error("nan score 0");
			}
			this->pSimulation()->score(pEffect->pEffectInfo(),
				this->pSimulation()->score(pEffect->pEffectInfo()) + score);
		}

		for (unsigned i = 0;
			i < this->pEndowmentFunction()->rEffects().size();
			i++)
		{
			Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];

			double score = 0;
			if (this->lpNetworkCache->outTieExists(alter))
			{
				score += this->lendowmentEffectContribution[alter][i];
			}

			if (curSetting) {
				permIter->reset();
			}

			// subtract sum over all permitted
			int j = 0;
			for (int alteri = 0; alteri < m; alteri++)
			{
				j = alteri;
				if (permIter)
				{
					if (!permIter->valid())
					{
						Rf_error("iterator not valid");
					}
					j = permIter->actor();
					permIter->next();
				}
				if (this->lpNetworkCache->outTieExists(j) && this->lpermitted[j])
				{
					score -= this->lendowmentEffectContribution[j][i]
						* this->lprobabilities[j];
				}
			}

			this->pSimulation()->score(pEffect->pEffectInfo(),
				this->pSimulation()->score(pEffect->pEffectInfo()) + score);
		}


		for (unsigned i = 0; i < this->pCreationFunction()->rEffects().size(); i++)
		{
			Effect * pEffect = this->pCreationFunction()->rEffects()[i];
			double score = 0;

			if (!this->lpNetworkCache->outTieExists(alter))
			{
				score += this->lcreationEffectContribution[alter][i];
			}

			if (curSetting) {
				permIter->reset();
			}

			// subtract sum over all permitted
			int j = 0;
			for (int alteri = 0; alteri < m; alteri++)
			{
				j = alteri;
				if (permIter)
				{
					if (!permIter->valid())
					{
						Rf_error("iterator not valid");
					}
					j = permIter->actor();
					permIter->next();
				}
				if (!this->lpNetworkCache->outTieExists(j) && this->lpermitted[j])
				{
					score -=
						this->lcreationEffectContribution[j][i] *
						this->lprobabilities[j];
				}
			}
			this->pSimulation()->score(pEffect->pEffectInfo(),
				this->pSimulation()->score(pEffect->pEffectInfo()) + score);
		}
	}
	if (permIter != 0)
	{
		delete permIter;
	}
}


// ----------------------------------------------------------------------------
// Section: symmetric networks methods
// ----------------------------------------------------------------------------

/**
 * Checks whether the alter would like to create the proposed link to the
 * current ego.
 */
void NetworkVariable::checkAlterAgreement(int alter)
{
	this->pSimulation()->pCache()->initialize(alter);
	this->preprocessEgo(alter);

	this->calculateSymmetricTieFlipContributions(this->lego, 1);

	this->calculateSymmetricTieFlipProbabilities(this->lego, 1, true);

	double probability = 0;
	double logprob = this->lsymmetricProbabilities[1];

	if (logprob > 0)
	{
		probability = 1.0 / (1.0 + exp(-logprob));
	}
	else
	{
		probability = exp(logprob);
		probability = probability / (1.0 + probability);
	}

	this->lsymmetricProbability = probability;
}

/**
 * Updates the scores for evaluation and endowment function effects according
 * to the current step in the simulation.
 */
void NetworkVariable::addAlterAgreementScores(bool accept)
{
	double probability = this->lsymmetricProbability;
	if (accept)
	{
		probability = 1 - probability;
	}

	for (unsigned i = 0;
		 i < this->pEvaluationFunction()->rEffects().size();
		 i++)
	{
		Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];

		double score = this->lsymmetricEvaluationEffectContribution[1][i]
			* probability;

		if (!accept)
		{
 			score = -1 * score;
 		}
		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}
	for (unsigned i = 0;
		 i < this->pEndowmentFunction()->rEffects().size();
		 i++)
	{
		Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];

		double score = 0;

		if (this->lpNetworkCache->outTieExists(this->lego))
		{
			score = this->lsymmetricEndowmentEffectContribution[1][i]
				* probability;
		}
		if (!accept)
 		{
 			score = -1 * score;
 		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}

	for (unsigned i = 0;
		i < this->pCreationFunction()->rEffects().size();
		i++)
	{
		Effect * pEffect = this->pCreationFunction()->rEffects()[i];

		double score = 0;

		if (!this->lpNetworkCache->outTieExists(this->lego))
		{
			score = this->lsymmetricCreationEffectContribution[1][i]
				* probability;
		}

		if (!accept)
 		{
 			score = -score;
 		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}
}
/**
 * Updates the scores for evaluation and endowment function effects according
 * to the current step in the simulation.
 */

void NetworkVariable::accumulateSymmetricModelScores(int alter, bool accept)
{
	double score = 0;
	double prEgo = 0;
	double prAlter = 0;
	double prSum = 0;

	switch(this->networkModelType())
	{
	case BFORCE:
		prEgo = this->lsymmetricProbabilities[0];
		for (unsigned i = 0;
			 i < this->pEvaluationFunction()->rEffects().size();
			 i++)
		{
			Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
			if (accept)
			{
				score = this->lsymmetricEvaluationEffectContribution[0][i]
					* (1 - prEgo);
			}
			else
			{
				score = -1 * this->lsymmetricEvaluationEffectContribution[0][i]
					* prEgo;
			}
			this->pSimulation()->score(pEffect->pEffectInfo(),
				this->pSimulation()->score(pEffect->pEffectInfo()) + score);
		}

		for (unsigned i = 0;
			 i < this->pEndowmentFunction()->rEffects().size();
			 i++)
		{
			Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];

			if (this->lpNetworkCache->outTieExists(alter))
			{
				if (accept)
				{
					score =
						this->lsymmetricEndowmentEffectContribution[0][i]
						* (1 - prEgo);
				}
				else
				{
					score =  -1  * prEgo *
						this->lsymmetricEndowmentEffectContribution[0][i];
				}
				this->pSimulation()->score(pEffect->pEffectInfo(),
					this->pSimulation()->score(pEffect->pEffectInfo()) + score);
			}
		}

		for (unsigned i = 0;
			 i < this->pCreationFunction()->rEffects().size();
			 i++)
		{
			Effect * pEffect = this->pCreationFunction()->rEffects()[i];

			if (!this->lpNetworkCache->outTieExists(alter))
			{
				if (accept)
				{
					score =
						this->lsymmetricCreationEffectContribution[0][i] *
						(1 - prEgo);
				}
				else
				{
					score =  -prEgo *
						this->lsymmetricCreationEffectContribution[0][i];
				}
				this->pSimulation()->score(pEffect->pEffectInfo(),
					this->pSimulation()->score(pEffect->pEffectInfo()) + score);
			}
		}

		break;

	case BAGREE:
		prEgo = this->lsymmetricProbabilities[0];
		prAlter = this->lsymmetricProbabilities[1];

		for (unsigned i = 0;
			 i < this->pEvaluationFunction()->rEffects().size();
			 i++)
		{
			Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
			if (!this->lpNetworkCache->outTieExists(alter))
			{
				if (accept)
				{
					score = this->lsymmetricEvaluationEffectContribution[0][i]
						* (1 - prEgo)  +
						this->lsymmetricEvaluationEffectContribution[1][i]
						* (1 - prAlter);
				}
				else
				{
					score =  -1  * prAlter * prEgo * ((1 - prEgo) *
						this->lsymmetricEvaluationEffectContribution[0][i] +
							(1 - prAlter) *
						this->lsymmetricEvaluationEffectContribution[1][i]) /
						(1 - prEgo * prAlter);
				}
			}
			else
			{
				if (accept)
				{
					score = (1 - prEgo) * (1 - prAlter) *
						(this->lsymmetricEvaluationEffectContribution[0][i]
							* prEgo +
							this->lsymmetricEvaluationEffectContribution[1][i]
							* prAlter) /
						(prEgo + prAlter - prEgo * prAlter);
				}
				else
				{
					score = -1 *
						((this->lsymmetricEvaluationEffectContribution[0][i]
							* prEgo) +
							(this->lsymmetricEvaluationEffectContribution[1][i]
								* prAlter));
				}
			}

			this->pSimulation()->score(pEffect->pEffectInfo(),
				this->pSimulation()->score(pEffect->pEffectInfo()) + score);
		}

		for (unsigned i = 0;
			 i < this->pEndowmentFunction()->rEffects().size();
			 i++)
		{
			Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];

			if (this->lpNetworkCache->outTieExists(alter))
			{
				if (accept)
				{
					score = (1 - prEgo) * (1 - prAlter) *
						(this->lsymmetricEndowmentEffectContribution[0][i]
							* prEgo + this->
							lsymmetricEndowmentEffectContribution[1][i]
							* prAlter) /
						(prEgo + prAlter - prEgo * prAlter);
				}
				else
				{
					score = -1 *
						((this->lsymmetricEndowmentEffectContribution[0][i]
							* prEgo) + (this->
								lsymmetricEndowmentEffectContribution[1][i]
								* prAlter));
				}

				this->pSimulation()->score(pEffect->pEffectInfo(),
					this->pSimulation()->score(pEffect->pEffectInfo()) + score);
			}
		}

		for (unsigned i = 0;
			 i < this->pCreationFunction()->rEffects().size();
			 i++)
		{
			Effect * pEffect = this->pCreationFunction()->rEffects()[i];

			if (!this->lpNetworkCache->outTieExists(alter))
			{
				if (accept)
				{
					score = (1 - prEgo) * (1 - prAlter) *
						(this->lsymmetricCreationEffectContribution[0][i] *
							prEgo +
							this->lsymmetricCreationEffectContribution[1][i] *
							prAlter) /
						(prEgo + prAlter - prEgo * prAlter);
				}
				else
				{
					score = -1 *
						(this->lsymmetricCreationEffectContribution[0][i] *
							prEgo +
							this->lsymmetricCreationEffectContribution[1][i] *
							prAlter);
				}
				this->pSimulation()->score(pEffect->pEffectInfo(),
					this->pSimulation()->score(pEffect->pEffectInfo()) + score);
			}

		}

		break;
	case BJOINT:

		prEgo = this->lsymmetricProbabilities[0];
		prAlter = this->lsymmetricProbabilities[1];
		prSum = prEgo + prAlter;
		if (prSum > 0)
		{
			prSum = 1.0 / (1.0 + exp(-prSum));
		}
		else
		{
			prSum = exp(prSum);
			prSum = prSum / (1.0 + prSum);
		}
		if (!accept)
		{
			prSum = 1.0 - prSum;
		}
		for (unsigned i = 0;
			 i < this->pEvaluationFunction()->rEffects().size();
			 i++)
		{
			Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
			score =  (1 - prSum) *
				(this->lsymmetricEvaluationEffectContribution[0][i] +
					this->lsymmetricEvaluationEffectContribution[1][i]);

			if (!accept)
			{
				score = -1 * score;
			}

			this->pSimulation()->score(pEffect->pEffectInfo(),
				this->pSimulation()->score(pEffect->pEffectInfo()) + score);
		}

		for (unsigned i = 0;
			 i < this->pEndowmentFunction()->rEffects().size();
			 i++)
		{
			Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];

			if (this->lpNetworkCache->outTieExists(alter))
			{
				score = (1.0 - prSum) *
					(this->lsymmetricEndowmentEffectContribution[0][i] +
						this->lsymmetricEndowmentEffectContribution[1][i]);

				if (!accept)
				{
					score = -1 * score;
				}
				this->pSimulation()->score(pEffect->pEffectInfo(),
					this->pSimulation()->score(pEffect->pEffectInfo()) + score);
			}
		}

		if (!this->lpNetworkCache->outTieExists(alter))
		{
			for (unsigned i = 0;
				 i < this->pCreationFunction()->rEffects().size();
				 i++)
			{
				Effect * pEffect = this->pCreationFunction()->rEffects()[i];

				score = (1 - prSum) *
					(this->lsymmetricCreationEffectContribution[0][i] +
						this->lsymmetricCreationEffectContribution[1][i]);

				if (!accept)
				{
					score = -score;
				}

				this->pSimulation()->score(pEffect->pEffectInfo(),
					this->pSimulation()->score(pEffect->pEffectInfo()) +
						score);
			}
		}

		break;

	case NORMAL:
	case AFORCE:
	case AAGREE:
	case DOUBLESTEP25:
	case DOUBLESTEP50:
	case DOUBLESTEP75:
	case DOUBLESTEP100:
	case NOTUSED:
		break;
	}
}

/**
 * For the given alter, this method calculates the contribution of each
 * effect if a tie from the ego to the alter
 * was flipped. These contributions are stored in arrays
 * <code>lsymmetricEvaluationEffectContribution</code>,
 * <code>lsymmetricEndowmentEffectContribution</code>, and
 * <code>lsymmetricCreationEffectContribution</code>. The first entry of the
 * array is for the ego effects, the second for the alter, controlled by the
 * integer parameter sub (0 or 1), as we need to do it both ways round.
 * Always permitted, never bipartite. Since always permitted we don't need
 * nan's to indicate not permitted.
 */
void NetworkVariable::calculateSymmetricTieFlipContributions(int alter,
int sub)
{
	int evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	int endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();
	int creationEffectCount = this->pCreationFunction()->rEffects().size();
	const vector<Effect *> & rEvaluationEffects =
		this->pEvaluationFunction()->rEffects();
	const vector<Effect *> & rEndowmentEffects =
		this->pEndowmentFunction()->rEffects();
	const vector<Effect *> & rCreationEffects =
		this->pCreationFunction()->rEffects();

	for (int i = 0; i < evaluationEffectCount; i++)
	{
		NetworkEffect * pEffect =
			(NetworkEffect *) rEvaluationEffects[i];
		double contribution =
			pEffect->calculateContribution(alter);

		// Tie withdrawals contribute in the opposite way

		if (this->lpNetworkCache->outTieExists(alter))
 		{
 			contribution = -contribution;
 		}
		this->lsymmetricEvaluationEffectContribution[sub][i] = contribution;
	}

	// The endowment effects have non-zero contributions on tie
	// withdrawals only. The opposite is true for tie creation effects.

	if (this->lpNetworkCache->outTieExists(alter) )
	{
		for (int i = 0; i < endowmentEffectCount; i++)
		{
			NetworkEffect * pEffect =
				(NetworkEffect *) rEndowmentEffects[i];
			this->lsymmetricEndowmentEffectContribution[sub][i] =
				-pEffect->calculateContribution(alter);
		}

		for (int i = 0; i < creationEffectCount; i++)
		{
			this->lsymmetricCreationEffectContribution[sub][i] = 0;
		}
	}
	else
	{
		for (int i = 0; i < creationEffectCount; i++)
		{
			NetworkEffect * pEffect =
				(NetworkEffect *) rCreationEffects[i];
			this->lsymmetricCreationEffectContribution[sub][i] =
				pEffect->calculateContribution(alter);
		}

		for (int i = 0; i < endowmentEffectCount; i++)
		{
			this->lsymmetricEndowmentEffectContribution[sub][i] = 0;
		}
	}

}

/**
 * Calculates the linear combinations which will be used to create the
 * probabilities of accepting change for symmetric networks.
 */
void NetworkVariable::calculateSymmetricTieFlipProbabilities(int alter,
													int sub, bool aagree)
{
	NetworkLongitudinalData * pData = (NetworkLongitudinalData *) this->pData();
	int evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	int endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();
	int creationEffectCount = this->pCreationFunction()->rEffects().size();

	// Calculate the total contribution of all effects

	double contribution = 0;

	for (int i = 0; i < evaluationEffectCount; i++)
	{
		Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
		contribution +=
			pEffect->parameter() *
			this->lsymmetricEvaluationEffectContribution[sub][i];
	}

	if (this->lpNetworkCache->outTieExists(alter))
	{
		for (int i = 0; i < endowmentEffectCount; i++)
		{
			Effect * pEffect =
				this->pEndowmentFunction()->rEffects()[i];
			contribution +=
				pEffect->parameter() *
				this->lsymmetricEndowmentEffectContribution[sub][i];
		}
	}
	else
	{

		for (int i = 0; i < creationEffectCount; i++)
		{
			Effect * pEffect = this->pCreationFunction()->rEffects()[i];
			contribution +=
				pEffect->parameter() *
					this->lsymmetricCreationEffectContribution[sub][i];
		}
	}

	if (aagree && (sub==1))
	{
		this->lsymmetricProbabilities[sub] = contribution +
										pData->universalOffset();
	}
	else
	{
		this->lsymmetricProbabilities[sub] = contribution;
	}
}
/**
 * Proposes and calculates probabilities for change.
 * Used for symmetric networks with model types beginning with B only.
 */
bool NetworkVariable::calculateModelTypeBProbabilities()
{
	this->preprocessEgo(this->lego);
	this->calculatePermissibleChanges();

	// choose alter
	int alter = this->lego;

	double * cumulativeRates = new double[this->n()];

	int numberPermitted = 0;
	for (int i = 0; i < this->n(); i++)
	{
		if (this->lpermitted[i] && i != this->lego)
		{
			numberPermitted++;
		}
		cumulativeRates[i] = this->rate(i);
		if (i > 0)
		{
			cumulativeRates[i] += cumulativeRates[i - 1];
		}

    }

	if (numberPermitted > 1)
	{
		while (alter == this->lego)
		{
			alter = nextIntWithCumulativeProbabilities(this->n(),
				cumulativeRates);
		}
	}

	this->lalterProbability = this->rate(alter) / cumulativeRates[this->n() - 1];

	delete [] cumulativeRates;

	this->lalter = alter;
//	Rprintf("net ego %d alter %d\n",this->lego, this->lalter);

	if (numberPermitted == 0 ||  (!this->lpermitted[alter]) ||
		alter == this->lego)
	{
		return false;
	}

	// calculate the probabilities

	double probability = 0;
	double sumValue = 0;
	this->pSimulation()->pCache()->initialize(alter);
	this->preprocessEgo(alter);
	this->calculateSymmetricTieFlipContributions(this->lego, 1);
	this->calculateSymmetricTieFlipProbabilities(this->lego, 1, false);

	this->pSimulation()->pCache()->initialize(this->lego);
	this->preprocessEgo(this->lego);
	this->calculateSymmetricTieFlipContributions(alter, 0);
	this->calculateSymmetricTieFlipProbabilities(alter, 0, false);
	switch(this->networkModelType())
	{
	case BFORCE:
		if (this->lsymmetricProbabilities[0] > 0)
		{
			probability = 1.0 / ( 1.0 + exp(-this->lsymmetricProbabilities[0]));
		}
		else
		{
			probability = exp(this->lsymmetricProbabilities[0]);
			probability = probability / (1.0 + probability);
		}
		this->lsymmetricProbabilities[0] = probability;
		break;

	case BAGREE:
		if (this->lsymmetricProbabilities[0] > 0)
		{
			probability = 1.0 / (1.0 + exp(-this->lsymmetricProbabilities[0]));
		}
		else
		{
			probability = exp(this->lsymmetricProbabilities[0]);
			probability = probability / (1.0 + probability);
		}
		this->lsymmetricProbabilities[0] = probability;

		if (this->lsymmetricProbabilities[1] > 0)
		{
			probability = 1.0 / (1.0 + exp(this->lsymmetricProbabilities[1]));

		}
		else
		{
			probability = exp(this->lsymmetricProbabilities[1]);
			probability = probability / (1.0 + probability);
		}
		this->lsymmetricProbabilities[1] = probability;

		if (!this->lpNetworkCache->outTieExists(alter))
		{
			probability = this->lsymmetricProbabilities[0] *
				this->lsymmetricProbabilities[1];
		}
		else
		{
			probability = this->lsymmetricProbabilities[0] +
				this->lsymmetricProbabilities[1] -
				this->lsymmetricProbabilities[0] *
				this->lsymmetricProbabilities[1];
		}
		break;

	case BJOINT:
		sumValue = this->lsymmetricProbabilities[0] +
			this->lsymmetricProbabilities[1];
		if (sumValue > 0)
		{
			probability = 1.0 / (1.0 + exp(-this->lsymmetricProbabilities[0] -
					this->lsymmetricProbabilities[1]));
		}
		else
		{
			probability = exp(this->lsymmetricProbabilities[0] +
				this->lsymmetricProbabilities[1]);
			probability = probability / (1.0 + probability);
		}
		break;

	case NORMAL:
	case AFORCE:
	case AAGREE:
	case NOTUSED:
	case DOUBLESTEP25:
	case DOUBLESTEP50:
	case DOUBLESTEP75:
	case DOUBLESTEP100:
		break;
	}

	this->lsymmetricProbability = probability;

	return true;
	// accept = nextDouble() < probability;

	// if (this->pSimulation()->pModel()->needScores())
	// {
	// 	this->accumulateSymmetricModelScores(alter, accept);
	// }

	// return accept;
}

// ----------------------------------------------------------------------------
// Section: Maximum likelihood related methods
// ----------------------------------------------------------------------------

/**
 * Calculates the probability of the given ministep assuming that the
 * ego of the ministep will change this variable.
 */
double NetworkVariable::probability(MiniStep * pMiniStep)
{
	// Initialize the cache object for the current ego
	this->pSimulation()->pCache()->initialize(pMiniStep->ego());

	NetworkChange * pNetworkChange =
		dynamic_cast<NetworkChange *>(pMiniStep);
	this->lego = pNetworkChange->ego();
	if (this->symmetric() && this->networkModelTypeB())
	{
		this->calculateModelTypeBProbabilities();
		if (this->pSimulation()->pModel()->needScores())
		{
			this->accumulateSymmetricModelScores(pNetworkChange->alter(),
				!pNetworkChange->diagonal());
		}
		// TODO no derivatives for symmetric models yet
		// if (this->pSimulation()->pModel()->needDerivatives())
		// {
		// 	this->accumulateDerivatives();
		// }
	}
	else
	{
		this->calculateTieFlipProbabilities();
		if (this->pSimulation()->pModel()->needScores())
		{
			this->accumulateScores(pNetworkChange->alter());
		}
		if (this->pSimulation()->pModel()->needDerivatives())
		{
			this->accumulateDerivatives();
		}
	}
	return this->lprobabilities[pNetworkChange->alter()];
}

/**
 * Updates the derivatives for effects
 * according to the current miniStep in the chain.
 */
void NetworkVariable::accumulateDerivatives() const
{
	int totalEvaluationEffects = this->pEvaluationFunction()->rEffects().size();
	int totalEndowmentEffects = this->pEndowmentFunction()->rEffects().size();
	int totalCreationEffects = this->pCreationFunction()->rEffects().size();
	int totalEffects =
		totalEvaluationEffects + totalEndowmentEffects + totalCreationEffects;
	Effect * pEffect1;
	Effect * pEffect2;
	double derivative;
	double * product = new double[totalEffects];
	double contribution1 = 0.0;
	double contribution2 = 0.0;

//	Rprintf("%d %d %d\n",totalEvaluationEffects,totalEndowmentEffects,
//		totalEffects);
	for (int effect1 = 0; effect1 < totalEffects; effect1++)
	{
		product[effect1] = 0.0;
		int endowment1 = effect1 - totalEvaluationEffects;
		int creation1 = effect1 - totalEvaluationEffects - totalEndowmentEffects;

		if (effect1 < totalEvaluationEffects)
		{
			pEffect1 = this->pEvaluationFunction()->rEffects()[effect1];
		}
		else if (effect1 < totalEvaluationEffects + totalEndowmentEffects)
		{
			pEffect1 = this->pEndowmentFunction()->rEffects()[endowment1];
		}
		else
		{
			pEffect1 = this->pCreationFunction()->rEffects()[creation1];
		}

		for (int alter = 0; alter < this->m(); alter++)
		{
			//	Rprintf("accum %d %d %d\n", alter, this->lpermitted[alter], this->lego);
			if (this->lpermitted[alter])
			{
				if (effect1 < totalEvaluationEffects)
				{
					product[effect1] +=
						this->levaluationEffectContribution[alter][effect1] *
						this->lprobabilities[alter];
				}
				else if (effect1 < totalEvaluationEffects + totalEndowmentEffects)
				{
					product[effect1] +=
						this->lendowmentEffectContribution[alter][endowment1] *
						this->lprobabilities[alter];
				}
				else
				{
					product[effect1] +=
						this->lcreationEffectContribution[alter][creation1] *
						this->lprobabilities[alter];
				}

				//	Rprintf("%d %d %f\n", alter, effect1, product[effect1]);
			}
		}
		for (int effect2 = effect1; effect2 < totalEffects; effect2++)
		{
			int endowment2 = effect2 - totalEvaluationEffects;
			int creation2 = effect2 - totalEvaluationEffects -
				totalEndowmentEffects;
			derivative = 0.0;

			if (effect2 < totalEvaluationEffects)
			{
				pEffect2 = this->pEvaluationFunction()->rEffects()[effect2];
			}
			else if (effect2 < totalEvaluationEffects + totalEndowmentEffects)
			{
				pEffect2 =
					this->pEndowmentFunction()->rEffects()[endowment2];
			}
			else
			{
				pEffect2 = this->pCreationFunction()->rEffects()[creation2];
			}

			for (int alter = 0; alter < this->m(); alter++)
			{
				if (this->lpermitted[alter])
				{
					if (effect1 < totalEvaluationEffects)
					{
						contribution1 =
							this->levaluationEffectContribution[alter][effect1];
					}
					else if (effect1 <
						totalEvaluationEffects + totalEndowmentEffects)
					{
						contribution1 =
							this->lendowmentEffectContribution[alter][endowment1];
					}
					else
					{
						contribution1 =
							this->lcreationEffectContribution[alter][creation1];
					}

					if (effect2 < totalEvaluationEffects)
					{
						contribution2 =
							this->levaluationEffectContribution[alter][effect2];
					}
					else if (effect2 <
						totalEvaluationEffects + totalEndowmentEffects)
					{
						contribution2 =
							this->lendowmentEffectContribution[alter][endowment2];
					}
					else
					{
						contribution2 =
							this->lcreationEffectContribution[alter][creation2];
					}


					derivative -=
						contribution1 * contribution2 *
						this->lprobabilities[alter];
					// Rprintf("deriv 2 %d %d %d %d %d %f %f %f %f %x %x\n",
					//alter,
					// 	effect1,
					// 	effect2, relativeEffect1, relativeEffect2,
					// 	derivative, contribution1, contribution2,
					// 	//	this->levaluationEffectContribution[alter][effect1],
					// 	//this->levaluationEffectContribution[alter][effect2],
					// 	this->lprobabilities[alter], pEffect1, pEffect2);
				}
			}

			this->pSimulation()->derivative(pEffect1->pEffectInfo(),
				pEffect2->pEffectInfo(),
				this->pSimulation()->derivative(pEffect1->pEffectInfo(),
					pEffect2->pEffectInfo()) +	derivative);
		}
	}

	for (int effect1 = 0; effect1 < totalEffects; effect1++)
	{
		int endowment1 = effect1 - totalEvaluationEffects;
		int creation1 = effect1 - totalEvaluationEffects - totalEndowmentEffects;

		for (int effect2 = effect1; effect2 < totalEffects; effect2++)
		{
			int endowment2 = effect2 - totalEvaluationEffects;
			int creation2 = effect2 - totalEvaluationEffects -
				totalEndowmentEffects;
			if (effect1 < totalEvaluationEffects)
			{
				pEffect1 = this->pEvaluationFunction()->rEffects()[effect1];
			}
			else if (effect1 < totalEvaluationEffects + totalEndowmentEffects)
			{
				pEffect1 = this->pEndowmentFunction()->rEffects()[endowment1];
			}
			else
			{
				pEffect1 = this->pCreationFunction()->rEffects()[creation1];
			}

			if (effect2 < totalEvaluationEffects)
			{
				pEffect2 = this->pEvaluationFunction()->rEffects()[effect2];
			}
			else if (effect2 < totalEvaluationEffects + totalEndowmentEffects)
			{
				pEffect2 = this->pEndowmentFunction()->rEffects()[endowment2];
			}
			else
			{
				pEffect2 = this->pCreationFunction()->rEffects()[creation2];
			}

			this->pSimulation()->derivative(pEffect1->pEffectInfo(),
				pEffect2->pEffectInfo(),
				this->pSimulation()->derivative(pEffect1->pEffectInfo(),
						pEffect2->pEffectInfo()) +
					product[effect1] * product[effect2]);
		}
	}

	delete [] product;
}

//	Rprintf("deriv %f\n", derivative;



/**
 * Returns whether applying the given ministep on the current state of this
 * variable would be valid with respect to all constraints. One can disable
 * the checking of up-only and down-only conditions.
 */
bool NetworkVariable::validMiniStep(const MiniStep * pMiniStep,
	bool checkUpOnlyDownOnlyConditions) const
{
	bool valid = DependentVariable::validMiniStep(pMiniStep);

	if (valid && !pMiniStep->diagonal())
	{
		NetworkLongitudinalData * pData =
			(NetworkLongitudinalData *) this->pData();
		const NetworkChange * pNetworkChange =
			dynamic_cast<const NetworkChange *>(pMiniStep);
		int i = pNetworkChange->ego();
		int j = pNetworkChange->alter();

		if (this->lpNetwork->tieValue(i, j))
		{
			if (checkUpOnlyDownOnlyConditions)
			{
				valid = !pData->upOnly(this->period());
			}
		}
		else
		{
			if (checkUpOnlyDownOnlyConditions)
			{
				valid = !pData->downOnly(this->period());
			}

			valid &= this->lpNetwork->outDegree(i) < pData->maxDegree() &&
				this->pReceivers()->active(j);
		}

		if (valid)
		{
			valid = !pData->structural(i, j, this->period());
		}

		// The filters may add some more conditions.

		for (unsigned i = 0;
			i < this->lpermittedChangeFilters.size() && valid;
			i++)
		{
			PermittedChangeFilter * pFilter = this->lpermittedChangeFilters[i];
			valid = pFilter->validMiniStep(pNetworkChange);
		}
	}

	return valid;
}


/**
 * Generates a random ministep for the given ego.
 */
MiniStep * NetworkVariable::randomMiniStep(int ego)
{
	this->pSimulation()->pCache()->initialize(ego);
	this->lego = ego;
	this->calculateTieFlipProbabilities();

	int m = 0;
	if (this->oneModeNetwork())
	{
		m = this->m();
	}
	else
	{
		m = this->m() + 1;
	}
	int alter = nextIntWithProbabilities(m, this->lprobabilities);

	MiniStep * pMiniStep =
		new NetworkChange(this->lpData, ego, alter,
			this->diagonalMiniStep(ego, alter));
	pMiniStep->logChoiceProbability(log(this->lprobabilities[alter]));

	return pMiniStep;
}


/**
 * Returns if the observed value for the option of the given ministep
 * is missing at either end of the period.
 */
bool NetworkVariable::missing(const MiniStep * pMiniStep) const
{
	const NetworkChange * pNetworkChange =
		dynamic_cast<const NetworkChange *>(pMiniStep);

	return this->lpData->missing(pNetworkChange->ego(),
		pNetworkChange->alter(),
		this->period()) ||
		this->lpData->missing(pNetworkChange->ego(),
			pNetworkChange->alter(),
			this->period() + 1);
}

/**
 * Returns if the given ministep is structurally determined in the period.
 */
bool NetworkVariable::structural(const MiniStep * pMiniStep) const
{
	const NetworkChange * pNetworkChange =
		dynamic_cast<const NetworkChange *>(pMiniStep);
	return !pMiniStep->diagonal() &&
		this->lpData->structural(pNetworkChange->ego(),
			pNetworkChange->alter(),
			this->period());
}


// ----------------------------------------------------------------------------
// Section: Properties
// ----------------------------------------------------------------------------

/**
 * Returns if this is a network variable.
 */
bool NetworkVariable::networkVariable() const
{
	return true;
}

/**
 * Returns if this is a symmetric network variable.
 */
bool NetworkVariable::symmetric() const
{
	const OneModeNetworkLongitudinalData * pData =
		dynamic_cast<const OneModeNetworkLongitudinalData *>(this->lpData);
	if (pData)
	{
		return pData->symmetric();
	}
	else
	{
		return false;
	}
}

/**
 * Returns if there are any constraints on the permitted changes of this
 * variable.
 */
bool NetworkVariable::constrained() const
{
	return DependentVariable::constrained() ||
		!this->lpermittedChangeFilters.empty();
}

/**
 * Returns the value of the alter in the current step. Only used for model type
 * B with symmetric networks; and for model type DOUBLESTEP
 */
int NetworkVariable::alter() const
{
//	Rprintf("in net %d\n", this->lalter);
	return this->lalter;
}

/**
 * Returns whether this is a diagonal step, for non-symmetric cases only
 */
bool NetworkVariable::diagonalMiniStep(int ego, int alter) const
{
	return (!this->oneModeNetwork() && alter == this->m()) ||
		(this->oneModeNetwork() && ego == alter);
}

void NetworkVariable::initializeSetting() {
	NetworkLongitudinalData * pNetworkData =
		dynamic_cast<NetworkLongitudinalData *>(this->pData());
	Setting* setting = lsettings[stepType()];
	string covariateName = pNetworkData->rSettingNames().at(this->stepType()).getCovarName();
	if (pSimulation()->pData()->pConstantDyadicCovariate(covariateName)) {
		setting->initDyadicSetting(
				pSimulation()->pData()->pConstantDyadicCovariate(covariateName)->rRowValues(
					ego()), ego());
	}
	if (pSimulation()->pData()->pChangingDyadicCovariate(covariateName)) {
		setting->initDyadicSetting(
				pSimulation()->pData()->pChangingDyadicCovariate(covariateName)->rRowValues(
					ego(), period()), ego());
	}
	setting->initSetting(ego());
}

const Setting * NetworkVariable::setting(int i) const {
	if (i >= numberSettings()) return 0;
	return lsettings[i];
}

}
