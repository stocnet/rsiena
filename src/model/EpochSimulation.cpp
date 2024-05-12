/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EpochSimulation.cpp
 *
 * Description: This file contains the implementation of the
 * EpochSimulation class.
 *****************************************************************************/

#include <algorithm>
#include <cmath>
#include <sstream>
#include <R_ext/Error.h>
#include <R_ext/Print.h>
#include <R_ext/Arith.h>
#include "SdeSimulation.h"
#include "EpochSimulation.h"
#include "utils/Random.h"
#include "utils/Utils.h"
#include "data/ActorSet.h"
#include "data/ExogenousEvent.h"
#include "data/LongitudinalData.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ContinuousLongitudinalData.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/ContinuousVariable.h"
#include "model/Model.h"
#include "model/SimulationActorSet.h"
#include "model/State.h"
#include "model/effects/Effect.h"
#include "model/tables/Cache.h"
#include "model/filters/AtLeastOneFilter.h"
#include "model/filters/DisjointFilter.h"
#include "model/filters/HigherFilter.h"
#include "model/filters/LowerFilter.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include <Rinternals.h>

using namespace std;

namespace siena {

SEXP getMiniStepList(const MiniStep& miniStep, int period);
SEXP getMiniStepDF(const MiniStep& miniStep);

// ----------------------------------------------------------------------------
// Section: Constructors and destructors
// ----------------------------------------------------------------------------

/**
 * Creates an epoch simulation object for the given observed data and an
 * actor-based model for that data.
 */
EpochSimulation::EpochSimulation(Data * pData, Model * pModel) {
	this->lpData = pData;
	this->lpModel = pModel;
	this->lpConditioningVariable = 0;

	// Create a cache object to be used to speed up effect calculations
	// during the simulation.

	this->lpCache = new Cache();

	// Create a wrapper for each actor set for simulation purposes,
	// and find the maximum number of actors in any actor set.

	int maxN = 0;

	for (unsigned i = 0; i < pData->rActorSets().size(); i++) {
		const ActorSet * pActorSet = pData->rActorSets()[i];
		SimulationActorSet * pSimulationActorSet = new SimulationActorSet(
				pActorSet);

		this->lsimulationActorSets.push_back(pSimulationActorSet);
		this->lactorSetMap[pActorSet] = pSimulationActorSet;
		maxN = std::max(maxN, pActorSet->n());
	}

	// Create the dependent variables from the observed data

	for (unsigned i = 0; i < pData->rDependentVariableData().size(); i++) {

		NetworkLongitudinalData * pNetworkData =
				dynamic_cast<NetworkLongitudinalData *>(pData->rDependentVariableData()[i]);
		BehaviorLongitudinalData * pBehaviorData =
				dynamic_cast<BehaviorLongitudinalData *>(pData->rDependentVariableData()[i]);
		ContinuousLongitudinalData * pContinuousData =
			dynamic_cast<ContinuousLongitudinalData *>(pData->rDependentVariableData()[i]);

		DependentVariable * pVariable = 0;
		ContinuousVariable * pContinuousVariable = 0;

		if (pNetworkData) 
		{
			pVariable = new NetworkVariable(pNetworkData, this);
			this->lvariables.push_back(pVariable);
			this->lvariableMap[pVariable->name()] = pVariable;
		} 
		else if (pBehaviorData) 
		{
			pVariable = new BehaviorVariable(pBehaviorData, this);
			this->lvariables.push_back(pVariable);
			this->lvariableMap[pVariable->name()] = pVariable;
		} 
		else if (pContinuousData) 
		{
			pContinuousVariable = new ContinuousVariable(pContinuousData, this);
			this->lcontinuousVariables.push_back(pContinuousVariable);
			this->lcontinuousVariableMap[pContinuousVariable->name()]
				= pContinuousVariable;
		} 
		else 
		{
			throw logic_error(
					"EpochSimulation: Network or behavior data expected.");
		}

		if (pModel->conditional()
				&& pModel->conditionalDependentVariable()
						== pVariable->name()) {
			this->lpConditioningVariable = pVariable;
		}
	}

	// Initialize the rate, evaluation, endowment, and creation
	// functions of all variables.

	for (unsigned i = 0; i < this->lvariables.size(); i++) {
		this->lvariables[i]->initializeRateFunction();
		this->lvariables[i]->initializeEvaluationFunction();
		this->lvariables[i]->initializeEndowmentFunction();
		this->lvariables[i]->initializeCreationFunction();
	}

	// Initialize all the effects for the continuous variables.
	
	for (unsigned i = 0; i < this->lcontinuousVariables.size(); i++) {
		this->lcontinuousVariables[i]->initializeFunction();
	}
	
	// Add network constraints to network variables.

	for (unsigned i = 0; i < pData->rNetworkConstraints().size(); i++) {
		const NetworkConstraint * pConstraint = pData->rNetworkConstraints()[i];
		NetworkVariable * pVariable1 =
				dynamic_cast<NetworkVariable *>(this->lvariableMap[pConstraint->networkName1()]);
		NetworkVariable * pVariable2 =
				dynamic_cast<NetworkVariable *>(this->lvariableMap[pConstraint->networkName2()]);

		if (!pVariable1) {
			throw logic_error(
					"Network variable " + pConstraint->networkName1()
							+ " expected.");
		}

		if (!pVariable2) {
			throw logic_error(
					"Network variable " + pConstraint->networkName2()
							+ " expected.");
		}

		if (pConstraint->type() == HIGHER) {
			pVariable1->addPermittedChangeFilter(
					new HigherFilter(pVariable1, pVariable2));
			pVariable2->addPermittedChangeFilter(
					new LowerFilter(pVariable2, pVariable1));
		} else if (pConstraint->type() == DISJOINT) {
			pVariable1->addPermittedChangeFilter(
					new DisjointFilter(pVariable1, pVariable2));
			pVariable2->addPermittedChangeFilter(
					new DisjointFilter(pVariable2, pVariable1));
		} else if (pConstraint->type() == AT_LEAST_ONE) {
			pVariable1->addPermittedChangeFilter(
					new AtLeastOneFilter(pVariable1, pVariable2));
			pVariable2->addPermittedChangeFilter(
					new AtLeastOneFilter(pVariable2, pVariable1));
		} else {
			std::ostringstream ss;
			ss << "Unexpected constraint type " << pConstraint->type()
					<< std::endl;
			throw logic_error(ss.str());
		}
	}

	// Allocate a helper array

	this->lcummulativeRates = new double[std::max(maxN,
			(int) this->lvariables.size())];
	this->ltargetChange = 0;

	// Create an SDE model for the evolution of the continuous variables
	
	if (this->lcontinuousVariables.size() > 0) {
		this->lpSdeSimulation = new SdeSimulation(this);
	} else {
		this->lpSdeSimulation = 0;
	}
	
	// Create a state object that will store the current values of all
	// dependent variables during the simulation.

	this->lpState = new State(this);
	this->lpChain = new Chain(pData);
}

/**
 * Deallocates this simulation object.
 */
EpochSimulation::~EpochSimulation() {
	delete[] this->lcummulativeRates;
	delete this->lpState;
	delete this->lpCache;
	delete this->lpChain;
	delete this->lpSdeSimulation;

	this->lcummulativeRates = 0;
	this->lpState = 0;
	this->lpCache = 0;
	this->lpSdeSimulation = 0;

	deallocateVector(this->lvariables);
	deallocateVector(this->lcontinuousVariables);
	deallocateVector(this->lsimulationActorSets);

	this->lvariableMap.clear();
	this->lcontinuousVariableMap.clear();
}

// ----------------------------------------------------------------------------
// Section: Model simulations
// ----------------------------------------------------------------------------

/**
 * Initializes the dependent variables as of the beginning of the specified
 * period.
 */
void EpochSimulation::initialize(int period) {
	this->lperiod = period;

	// Initialize the active actor indicators

	for (unsigned i = 0; i < this->lsimulationActorSets.size(); i++) {
		SimulationActorSet * pActorSet = this->lsimulationActorSets[i];

		for (int i = 0; i < pActorSet->n(); i++) {
			pActorSet->active(i,
					this->lpData->active(pActorSet->pOriginalActorSet(), i,
							period));
		}
	}

	// Initialize each discrete dependent variable
	for (unsigned i = 0; i < this->lvariables.size(); i++) {
		this->lvariables[i]->initialize(period);
	}

	// Initialize each continuous dependent variable
	for (unsigned i = 0; i < this->lcontinuousVariables.size(); i++) {
		this->lcontinuousVariables[i]->initialize(period);
	}

	// Initialize the effects for the upcoming simulation

	for (unsigned i = 0; i < this->lvariables.size(); i++) {
		const Function * pFunction = this->lvariables[i]->pEvaluationFunction();

		for (unsigned j = 0; j < pFunction->rEffects().size(); j++) {
			pFunction->rEffects()[j]->initialize(this->lpData, this->lpState,
					period, this->lpCache);
		}

		pFunction = this->lvariables[i]->pEndowmentFunction();

		for (unsigned j = 0; j < pFunction->rEffects().size(); j++) {
			pFunction->rEffects()[j]->initialize(this->lpData, this->lpState,
					period, this->lpCache);
		}

		pFunction = this->lvariables[i]->pCreationFunction();

		for (unsigned j = 0; j < pFunction->rEffects().size(); j++) {
			pFunction->rEffects()[j]->initialize(this->lpData, this->lpState,
					period, this->lpCache);
		}
	}

	for (unsigned i = 0; i < this->lcontinuousVariables.size(); i++) {
		const Function * pFunction = this->lcontinuousVariables[i]->pFunction();

		for (unsigned j = 0; j < pFunction->rEffects().size(); j++) {
			pFunction->rEffects()[j]->initialize(this->lpData, this->lpState,
					period, this->lpCache);
		}
	}
	
	// Initialize the SDE for the upcoming simulation
	if (this->lcontinuousVariables.size() > 0) {
		 this->lpSdeSimulation->initialize(period);
	}

	// Reset the time
	this->ltime = 0;

	// Exogenous events
	this->lpEvents = this->lpData->pEventSet(period);
	this->lnextEvent = this->lpEvents->begin();

	// targets for conditional simulation
	if (this->lpModel->conditional()) {
		this->ltargetChange = this->lpModel->targetChange(this->lpData, period);
	} else {
		this->ltargetChange = 0;
	}

	// Reset scores
	this->lscores.clear();
	// Reset derivatives
	this->lderivatives.clear();
	// reset chain
	//this->lpChain->clear(); make the user do this
	this->lpChain->period(period);
}

/**
 * Simulates one complete period for the model and data
 */
void EpochSimulation::runEpoch(int period) {
	this->initialize(period);
	for (int nIter = 0;; nIter++) {
		this->runStep();

		if (this->lpModel->conditional()) {
			if (this->lpConditioningVariable->simulatedDistance()
					>= this->ltargetChange) {
				break;
			} else if (nIter > 1000000) {
#ifdef STANDALONE
				exit(1);
#endif
#ifndef STANDALONE
				Rf_error("%s %s", "Unlikely to terminate this epoch:",
						" more than 1000000 steps");
#endif
			}
		} else {
			if (this->ltime >= 1) {
				break;
			} else if (nIter > 1000000) {
#ifdef STANDALONE
				exit(1);
#endif
#ifndef STANDALONE
				Rf_error("%s %s", "Unlikely to terminate this epoch:",
						" more than 1000000 steps");
#endif
			}
		}
	}
	if (this->lpEvents->size()) {
		this->setLeaversBack();
	}
	if (this->pModel()->needChain()) {
		this->calculateRates();
		this->pChain()->finalReciprocalRate(1 / this->grandTotalRate());
	}
}

/**
 * Simulates a single step of the actor-oriented model.
 */
void EpochSimulation::runStep() {
	this->calculateRates();
	this->drawTimeIncrement();
	double nextTime = this->ltime + this->ltau;
	DependentVariable * pSelectedVariable = 0;
	int selectedActor = 0;

	if (this->lpModel->conditional() || nextTime < 1) 
	{
		if (this->reachedCompositionChange()) {
			this->makeNextCompositionChange();
			if (this->pModel()->needScores()) {
				// not done if parallel running: bug in Siena3
				if (!this->lpModel->parallelRun()) {
					this->accumulateRateScores(this->ltau);
				}
			}
		} else {
			this->ltime = nextTime;
			// SDE step 
			if (this->lcontinuousVariables.size() > 0) {
				this->lpSdeSimulation->setBergstromCoefficients(this->ltau);
				this->updateContinuousVariablesAndScores();
			}
			
			pSelectedVariable = this->chooseVariable();
			selectedActor = this->chooseActor(pSelectedVariable);

			this->lpCache->initialize(selectedActor);

			pSelectedVariable->makeChange(selectedActor);			
						
			if (pSelectedVariable->networkModelTypeDoubleStep())  // perhaps make a double step
			{ 
				double value = nextDouble();
				if (value < pSelectedVariable->networkDoubleStepProb())
				{
					int chosenAlter = pSelectedVariable->alter();
					if (chosenAlter != selectedActor) // type DOUBLESTEP implies one-mode
					{
						this->lpCache->initialize(chosenAlter);
						pSelectedVariable->makeChange(chosenAlter); // here the double step is taken
					}
				}
			}

			if (pSelectedVariable->successfulChange()) {
				if (this->pModel()->needChain()) {
					// update rate probabilities on the current final ministep
					this->lpChain->pLast()->pPrevious()->logOptionSetProbability(
							log(pSelectedVariable->rate(selectedActor)
											/ this->grandTotalRate()));
					this->lpChain->pLast()->pPrevious()->reciprocalRate(
							1.0 / this->grandTotalRate());
				}
			}
			// Update the scores for rate parameters
			if (this->pModel()->needScores()) {
				this->accumulateRateScores(this->ltau, pSelectedVariable,
						selectedActor);
			}
		}
	} else 
	{
		// Make sure we stop at 1.0 precisely.

		this->ltau = 1 - this->ltime;
		this->ltime = 1;

		// SDE step 
		if (this->lcontinuousVariables.size() > 0) {
			this->lpSdeSimulation->setBergstromCoefficients(this->ltau);
			this->updateContinuousVariablesAndScores();
		}
			
		// Update rate scores
		if (this->pModel()->needScores()) 
		{
			this->accumulateRateScores(this->ltau);
		}
	}
}

/**
 * Calculates the rates of change of each actor for each dependent variable and
 * the total rates of change for each variable summed over all actors.
 */
void EpochSimulation::calculateRates() {
	this->lgrandTotalRate = 0;

	for (unsigned i = 0; i < this->lvariables.size(); i++) {
		this->lvariables[i]->calculateRates();
		this->lgrandTotalRate += this->lvariables[i]->totalRate();
	}
}

/**
 * Generates an exponential variate tau with the sum ot total rates over all
 * dependent variables as the distribution parameter. It is used later to
 * increment the current time of the simulation.
 */
void EpochSimulation::drawTimeIncrement() {
	// use QAD if parallel running as other one uses 2 random numbers
	// also use QAD if STANDALONE (SienaProfile.cpp)
	double tau;

#ifndef STANDALONE
	if (this->lpModel->parallelRun()) {
#endif
		tau = nextExponentialQAD(this->lgrandTotalRate);
#ifndef STANDALONE
	} else {
		tau = nextExponential(this->lgrandTotalRate);
	}
#endif

	this->ltau = tau;
}

/**
 * Returns if the simulation has reached the time point of the next
 * exogenous event of composition change.
 */
bool EpochSimulation::reachedCompositionChange() const {
	return this->lnextEvent != this->lpEvents->end()
			&& (*this->lnextEvent)->time() <= this->ltime + this->ltau;
}

/**
 * Makes the current composition change and resets the time of this simulation
 * to the time of the composition change.
 */
void EpochSimulation::makeNextCompositionChange() {
	ExogenousEvent * pEvent = *this->lnextEvent;
	this->lnextEvent++;

	SimulationActorSet * pActorSet = this->lactorSetMap[pEvent->pActorSet()];

	if (pEvent->type() == JOINING) {
		pActorSet->active(pEvent->actor(), true);

		for (unsigned i = 0; i < this->lvariables.size(); i++) {
			this->lvariables[i]->actOnJoiner(pActorSet, pEvent->actor());
		}
	} else if (pEvent->type() == LEAVING) {
		pActorSet->active(pEvent->actor(), false);

		for (unsigned i = 0; i < this->lvariables.size(); i++) {
			this->lvariables[i]->actOnLeaver(pActorSet, pEvent->actor());
		}
	}

	this->ltau = pEvent->time() - this->ltime;
	this->ltime = pEvent->time();
}

/**
 * Resets the values for any actors who left the system during the current
 * period to their value at the start of the period. It will then not affect
 * the calculation of statistics. In fact resets values for all non active
 * actors.
 */
void EpochSimulation::setLeaversBack() {
	for (unsigned i = 0; i < this->lvariables.size(); i++) {
// 		for (EventSet::iterator  iter = this->lpEvents->begin();
// 			 iter!=this->lpEvents->end();
// 			 iter++)
// 		{
// 			ExogenousEvent * pEvent = *iter;

// 			if (pEvent->type() == LEAVING)
// 			{
// 				this->lvariables[i]->setLeaverBack(pEvent->pActorSet(),
// 					pEvent->actor());
// 			}
// 		}

		DependentVariable *pVariable = this->lvariables[i];
		const SimulationActorSet *pActorSet = pVariable->pActorSet();

		for (int j = 0; j < pVariable->n(); j++) {
			if (!pActorSet->active(j)) {
				pVariable->setLeaverBack(pActorSet, j);
			}
		}
	}
}

/**
 * Chooses one of the dependent varaibles randomly with probabilities
 * proportional to the total rate of each variable.
 */
DependentVariable * EpochSimulation::chooseVariable() const {
	int index = 0;

	if (this->lvariables.size() > 1) {
		for (unsigned i = 0; i < this->lvariables.size(); i++) {
			this->lcummulativeRates[i] = this->lvariables[i]->totalRate();

			if (i > 0) {
				this->lcummulativeRates[i] += this->lcummulativeRates[i - 1];
			}
		}

		index = nextIntWithCumulativeProbabilities(this->lvariables.size(),
				this->lcummulativeRates);
		//	Rprintf(" %d %f %f %f\n", index, this->lcummulativeRates[0],
		//this->lcummulativeRates[1],
		//  this->lcummulativeRates[2]);
	}

	return this->lvariables[index];
}

/**
 * Chooses a random actor with probabilities proportional to the rate of change
 * for the given variable.
 */
int EpochSimulation::chooseActor(const DependentVariable * pVariable) const {
	for (int i = 0; i < pVariable->n(); i++) {
		this->lcummulativeRates[i] = pVariable->rate(i);

		if (i > 0) {
			this->lcummulativeRates[i] += this->lcummulativeRates[i - 1];
		}
	}

	return nextIntWithCumulativeProbabilities(pVariable->n(),
			this->lcummulativeRates);
}

/**
 * Accumulates the scores for the rate parameters.
 */
void EpochSimulation::accumulateRateScores(double tau,
		const DependentVariable * pSelectedVariable, int selectedActor) {
	for (unsigned i = 0; i < this->lvariables.size(); i++) {

		if (this->lvariables[i]->symmetric() && this->lvariables[i]->networkModelTypeB()) {
			//	Rprintf("1got here %d\n",this->lvariables[i]->alter());
			this->lvariables[i]->accumulateRateScores(tau, pSelectedVariable,
					selectedActor, this->lvariables[i]->alter());
			//Rprintf("got here2\n");
		} else {
			//	Rprintf("got else %d\n",this->lvariables[i]->alter());
			this->lvariables[i]->accumulateRateScores(tau, pSelectedVariable,
					selectedActor);
		}
	}
}

/**
 * Updates the continuous variables based on the SDE and the scores of
 * the effects for these variables (in case the score function method
 * is used).
 */
void EpochSimulation::updateContinuousVariablesAndScores() {
	// up to now only for one continuous variable
	// function is never called if lcontVar's.size() == 0
	if (this->lcontinuousVariables.size() > 1) {
		Rf_error("EpochSimulation: Not more than one continuous variable.");
	}

	ContinuousVariable * pVariable = this->lcontinuousVariables[0];
	pVariable->calculateEffectContribution();
	
	vector<double> actorMeans(pVariable->n());
	vector<double> actorErrors(pVariable->n());
	
	for (int actor = 0; actor < pVariable->n(); actor++) 
	{
		// new value = deterministic part + random part
		actorMeans[actor] = pVariable->totalFunctionContribution(actor);
		actorErrors[actor] = this->lpSdeSimulation->randomComponent();
		double newValue = actorMeans[actor] + actorErrors[actor];
		pVariable->value(actor, newValue);
		// Rprintf("t = %f: actor %d has value %f\n", this->ltime, actor, newValue);
	}

	if (this->pModel()->needScores())
	{
		pVariable->accumulateScores(actorMeans, actorErrors, this->ltau);
	}
}

// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the data object underlying this simulation.
 */
const Data * EpochSimulation::pData() const {
	return this->lpData;
}

/**
 * Returns the actor-based model simulated by this simulation object.
 */
const Model * EpochSimulation::pModel() const {
	return this->lpModel;
}

/**
 * Returns the object responsible for the simulation of the continuous 
 * variables.
 */
const SdeSimulation * EpochSimulation::pSdeSimulation() const {
	return this->lpSdeSimulation;
}

/**
 * Returns the chain representing the events simulated by this simulation object.
 */
Chain * EpochSimulation::pChain() {
	return this->lpChain;
}

/**
 * Sets the chain representing the events simulated by this object to the
 * given chain.
 */
void EpochSimulation::pChain(Chain * pChain) {
	delete this->lpChain;
	this->lpChain = pChain;
}

/**
 * Clears the chain representing the events simulated by this object to the
 * given chain.
 */
void EpochSimulation::clearChain() {
	this->lpChain->clear();
}

/**
 * Returns the dependent variable with the given name if it exists;
 * otherwise 0 is returned.
 */
const DependentVariable * EpochSimulation::pVariable(string name) const {
	map<string, DependentVariable *>::const_iterator iter =
			this->lvariableMap.find(name);
	const DependentVariable * pVariable = 0;

	if (iter != this->lvariableMap.end()) {
		pVariable = iter->second;
	}

	return pVariable;
}

/**
 * Returns the dependent variable with the given name if it exists;
 * otherwise 0 is returned.
 */
const ContinuousVariable * EpochSimulation::pContinuousVariable(string name) const {
	map<string, ContinuousVariable *>::const_iterator iter =
			this->lcontinuousVariableMap.find(name);
	const ContinuousVariable * pVariable = 0;

	if (iter != this->lcontinuousVariableMap.end()) {
		pVariable = iter->second;
	}

	return pVariable;
}

/**
 * Returns a reference to the vector of discrete dependent variables.
 */
const vector<DependentVariable *> & EpochSimulation::rVariables() const {
	return this->lvariables;
}

/**
 * Returns a reference to the vector of continuous dependent variables.
 */
const vector<ContinuousVariable *> & EpochSimulation::rContinuousVariables() const {
	return this->lcontinuousVariables;
}


/**
 * Returns the wrapper actor set corresponding to the given original actor set.
 */
const SimulationActorSet * EpochSimulation::pSimulationActorSet(
		const ActorSet * pOriginalActorSet) const {
	map<const ActorSet *, SimulationActorSet *>::const_iterator iter =
			this->lactorSetMap.find(pOriginalActorSet);

	const SimulationActorSet * pSimulationActorSet = 0;

	if (iter != this->lactorSetMap.end()) {
		pSimulationActorSet = iter->second;
	}

	return pSimulationActorSet;
}

/**
 * Returns the currently simulated period.
 */
int EpochSimulation::period() const {
	return this->lperiod;
}

/**
 * Returns the time taken in the simulation.
 */
double EpochSimulation::time() const {
	return this->ltime;
}

/**
 * Returns the current score for the given effect. The scores are updated
 * in each ministep of the simulation.
 */
double EpochSimulation::score(const EffectInfo * pEffect) const {
	map<const EffectInfo *, double>::const_iterator iter = this->lscores.find(
			pEffect);
	double score = 0;

	if (iter != this->lscores.end()) {
		score = iter->second;
	}

	return score;
}

/**
 * Sets the score for the given effect to the given value.
 */
void EpochSimulation::score(const EffectInfo * pEffect, double value) {
	this->lscores[pEffect] = value;
}

/**
 * Sets the score of the scale parameter to the given value. 
 */
void EpochSimulation::basicScaleScore(double score)
{
	this->lpSdeSimulation->basicScaleScore(score);
}

/**
 * Returns the current derivatives for the given pair of effects.
 * The derivatives are updated for each ministep of a chain.
 */
double EpochSimulation::derivative(const EffectInfo * pEffect1,
		const EffectInfo * pEffect2) const {
	map<const EffectInfo *, map<const EffectInfo *, double> >::const_iterator iter =
			this->lderivatives.find(pEffect1);
	double derivative = 0;

	if (iter != this->lderivatives.end()) {
		const map<const EffectInfo *, double> effect1Map = iter->second;
		map<const EffectInfo *, double>::const_iterator iter2 = effect1Map.find(
				pEffect2);
		if (iter2 != effect1Map.end()) {
			derivative = iter2->second;
		}
	}

	return derivative;
}

/**
 * Calculates the likelihood from the chain ignoring the log factorial for
 * constant rates.
 */
double EpochSimulation::calculateLikelihood() const {
	//Rprintf("here\n");
	double sumLogOptionSetProbabilities = 0;
	double sumLogChoiceProbabilities = 0;
	double loglik = 0;
//	Rprintf(" %d\n", this->lpChain->ministepCount());
	// set up array to store counts of structurally active ministeps by variable
	int *counts = new int[this->lvariables.size()];
	for (unsigned i = 0; i < this->lvariables.size(); i++) {
		counts[i] = 0;
	}

	// create array to store number of actors by variable
	int *nActors = new int[this->lvariables.size()];
	for (unsigned i = 0; i < this->lvariables.size(); i++) {
		nActors[i] = this->lvariables[i]->n();
	}

	MiniStep *pMiniStep = this->lpChain->pFirst()->pNext();

//Rprintf("%d %x %x\n", this->lpChain->ministepCount(),
// pMiniStep, this->lpChain->pLast());

	while (pMiniStep != this->lpChain->pLast()) {
		DependentVariable * pVariable =
				this->lvariables[pMiniStep->variableId()];
		sumLogOptionSetProbabilities += pMiniStep->logOptionSetProbability();
		sumLogChoiceProbabilities += pMiniStep->logChoiceProbability();
		//	if (!R_finite(pMiniStep->logChoiceProbability()))
		//{
		//		Rf_PrintValue(getMiniStepDF(*pMiniStep));
		//	Rprintf(" epcoh %f %f\n", pMiniStep->logOptionSetProbability(),
		//	pMiniStep->logChoiceProbability() );
		//}
		if (!pVariable->structural(pMiniStep)) {
			counts[pMiniStep->variableId()]++;
		}

		pMiniStep = pMiniStep->pNext();
	}

	loglik += sumLogChoiceProbabilities;
	if (!R_finite(pMiniStep->logChoiceProbability())) {
		Rprintf("sum choice %f", loglik);
	}
	if (this->lsimpleRates) {
		for (unsigned i = 0; i < this->lvariables.size(); i++) {
			DependentVariable * pVariable = this->lvariables[i];
			double lambda = pVariable->basicRate();
			loglik += counts[i] * log(lambda) - lambda * pVariable->n()
					- this->lnFactorial(counts[i]);
			// if (!R_finite(loglik))
			// {
			// 	Rprintf("basic rate %f count %d log %f %f\n",lambda, counts[0],
			// 		log(lambda), loglik);
			// }
		}
	} else {
		loglik += sumLogOptionSetProbabilities
				+ normalDensity(1, this->lpChain->mu(),
						sqrt(std::max(0.00, this->lpChain->sigma2())), 1)
				+ log(this->lpChain->finalReciprocalRate());
		// if (!R_finite(loglik))
		// {
		// 	Rprintf("mu %f sigma2 %f final %f sumop %f loglik %f\n",
		// 		this->lpChain->mu(), this->lpChain->sigma2(),
		// 		this->lpChain->finalReciprocalRate(),
		// 		sumLogOptionSetProbabilities, loglik);
		// }
	}

	delete[] counts;
	delete[] nActors;
	return loglik;

}
/**
 * Calculates log factorial of its integer argument
 */
double EpochSimulation::lnFactorial(int a) const {
	int y;
	double z;
	if (a == 1)
		return 0;
	else {
		z = 0;

		for (y = 2; y <= a; y++)

			z = log(static_cast<double>(y)) + z;

		return z;
	}
}
/**
 * Updates parameters. Probably legacy code and not used.
 */
void EpochSimulation::updateParameters(int period) {
	Rprintf("ever used?\n");
	for (unsigned i = 0; i < this->lvariables.size(); i++) {
		this->lvariables[i]->updateBasicRate(period);
		this->lvariables[i]->updateEffectParameters();
	}

}
/**
 * Returns the current map of derivatives for the given effect.
 * The derivatives are updated for each ministep of a chain.
 */
map<const EffectInfo*, double> EpochSimulation::derivative(
		const EffectInfo * pEffect) const {
	map<const EffectInfo *, map<const EffectInfo *, double> >::const_iterator iter =
			this->lderivatives.find(pEffect);

	map<const EffectInfo *, double> effectMap;
	if (iter != this->lderivatives.end()) {
		const map<const EffectInfo *, double> effectMap = iter->second;
	}

	return effectMap;
}

/**
 * Sets the derivative for the given effects to the given value.
 */
void EpochSimulation::derivative(const EffectInfo * pEffect1,
		const EffectInfo * pEffect2, double value) {
	this->lderivatives[pEffect1][pEffect2] = value;
}

/**
 * Returns the cache object used to speed up the simulations.
 */
Cache * EpochSimulation::pCache() const {
	return this->lpCache;
}

/**
 * Returns the total rate over all dependent variables.
 */
double EpochSimulation::grandTotalRate() const {
	return this->lgrandTotalRate;
}
/**
 * Stores if simple rates should be used in simulations.
 */
void EpochSimulation::simpleRates(bool flag) {
	this->lsimpleRates = flag;
}

/**
 * Returns if simple rates should be used in simulations.
 */
bool EpochSimulation::simpleRates() const {
	return this->lsimpleRates;
}

}
