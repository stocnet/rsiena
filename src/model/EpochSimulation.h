/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EpochSimulation.h
 *
 * Description: This file contains the definition of the
 * EpochSimulation class.
 *****************************************************************************/

#ifndef EPOCHSIMULATION_H_
#define EPOCHSIMULATION_H_

#include <vector>
#include <map>
#include <string>
#include "data/Data.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class DependentVariable;
class ContinuousVariable;
class Model;
class ActorSet;
class EffectInfo;
class SimulationActorSet;
class State;
class Cache;
class Chain;
class MiniStep;
class SdeSimulation;


// ----------------------------------------------------------------------------
// Section: EpochSimulation class
// ----------------------------------------------------------------------------

/**
 * This class provides the functionality necessary for simulating a model
 * between two observations.
 */
class EpochSimulation
{
public:
	EpochSimulation(Data * pData, Model * pModel);
	virtual ~EpochSimulation();

	void initialize(int period);

	// Method of moments related
	void runEpoch(int period);

	// Accessors
	const Data * pData() const;
	const Model * pModel() const;
	const SdeSimulation * pSdeSimulation() const;
	const DependentVariable * pVariable(std::string name) const;
	const ContinuousVariable * pContinuousVariable(std::string name) const;
	const std::vector<DependentVariable *> & rVariables() const;
	const std::vector<ContinuousVariable *> & rContinuousVariables() const;
	const SimulationActorSet * pSimulationActorSet(
			const ActorSet * pOriginalActorSet) const;
	int period() const;
	double time() const;
	Cache * pCache() const;

	// Scores

	double score(const EffectInfo * pEffect) const;
	void score(const EffectInfo * pEffect, double value);
	void basicScaleScore(double score);
	std::map<const EffectInfo *, double>
		derivative(const EffectInfo * pEffect1) const;
	double derivative(const EffectInfo * pEffect1,
			const EffectInfo * pEffect2) const;
	void derivative(const EffectInfo * pEffect1, const EffectInfo * pEffect2,
			double value);
	Chain * pChain();
	void  pChain(Chain * pChain);
	void clearChain();
	void updateParameters(int period);
	double calculateLikelihood() const;

	void simpleRates(bool flag);
	bool simpleRates() const;

	double lnFactorial(int a) const;
protected:
	void calculateRates();
	double grandTotalRate() const;
	DependentVariable * chooseVariable() const;
	int chooseActor(const DependentVariable * pVariable) const;

    // A vector of discrete dependent variables with their current values
	std::vector<DependentVariable *> lvariables;

	// A vector of continuous dependent variables with their current values
    std::vector<ContinuousVariable *> lcontinuousVariables;

private:
	void runStep();
	void drawTimeIncrement();
	bool reachedCompositionChange() const;
	void makeNextCompositionChange();
	void setLeaversBack();
	void accumulateRateScores(double tau,
			const DependentVariable * pSelectedVariable = 0,
			int selectedActor = 0);
	void updateContinuousVariablesAndScores();

	// The observed data the model is based on
	Data * lpData;

	// The actor-based model to be simulated
	Model * lpModel;

	// The addition to the lpModel for the simulation of
	// continuous variables
	SdeSimulation * lpSdeSimulation;

	// A wrapper object per actor set for simulation purposes
	std::vector<SimulationActorSet *> lsimulationActorSets;

	// Stores the wrappers of each original actor set
	std::map<const ActorSet *, SimulationActorSet *> lactorSetMap;

	// The dependent variable for look-ups by variable names
	std::map<std::string, DependentVariable *> lvariableMap;

	// The continuous dependent variable for look-ups by variable names
    std::map<std::string, ContinuousVariable *> lcontinuousVariableMap;
	
	// The current period to be simulated
	int lperiod {};

	// An array of cummulative rates used for the random selection of
	// the dependent variable to change and the actor to make the change.

	double * lcummulativeRates {};

	// The total rate over all dependent variables
	double lgrandTotalRate {};

	// The current time of the simulation
	double ltime {};

	// The current increment of time of the simulation
	double ltau {};

	// A sorted set of exogenous events of composition change
	const EventSet * lpEvents;

	// An iterator to the next event still to be processed.
	EventSet::const_iterator lnextEvent;

	// Target amount of change for this period if we are using conditional simulation
	int ltargetChange {};

	// The dependent variable the simulation is conditioned upon
	DependentVariable * lpConditioningVariable;

	// Values of scores in this simulation: one for each selected effect,
	// including the rate effects, but excluding the basic rate effect.

	std::map<const EffectInfo *, double> lscores;
	std::map<const EffectInfo *, std::map <const EffectInfo *, double> > lderivatives;

	State * lpState;
	Cache * lpCache;

	Chain * lpChain;
	bool lsimpleRates {};

};

}

#endif /*EPOCHSIMULATION_H_*/
