/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkVariable.h
 *
 * Description: This file contains the definition of the
 * NetworkVariable class.
 *****************************************************************************/

#ifndef NETWORKVARIABLE_H_
#define NETWORKVARIABLE_H_

#include <vector>

#include "DependentVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class NetworkLongitudinalData;
class NetworkCache;
class PermittedChangeFilter;
class ITieIterator;

// ----------------------------------------------------------------------------
// Section: NetworkVariable class
// ----------------------------------------------------------------------------

/**
 * This class represents the state of a one-mode network variable.
 * @see DependentVariable
 */
class NetworkVariable : public DependentVariable
{
public:
	NetworkVariable(NetworkLongitudinalData * pData,
		EpochSimulation * pSimulation);
	virtual ~NetworkVariable();

	void addPermittedChangeFilter(PermittedChangeFilter * pFilter);

	const SimulationActorSet * pSenders() const;
	const SimulationActorSet * pReceivers() const;
	virtual int m() const;
	virtual LongitudinalData * pData() const;
	bool oneModeNetwork() const;
	virtual bool networkVariable() const;
	virtual bool constrained() const;
	virtual bool symmetric() const;

	void networkModelType(int type);
	virtual NetworkModelType networkModelType() const;
	virtual bool networkModelTypeB() const;
	virtual bool networkModelTypeDoubleStep() const;
	virtual double networkDoubleStepProb() const;

	virtual void initialize(int period);
	virtual bool canMakeChange(int actor) const;
	virtual void makeChange(int actor);
	virtual void actOnJoiner(const SimulationActorSet * pActorSet, int actor);
	virtual void actOnLeaver(const SimulationActorSet * pActorSet, int actor);
	virtual void setLeaverBack(const SimulationActorSet * pActorSet,
		int actor);
	Network * pNetwork() const;

	int ego() const;
	virtual int alter() const;

	virtual double probability(MiniStep * pMiniStep);
	virtual bool validMiniStep(const MiniStep * pMiniStep,
		bool checkUpOnlyDownOnlyConditions = true) const;
	virtual MiniStep * randomMiniStep(int ego);
	virtual bool missing(const MiniStep * pMiniStep) const;
	virtual bool structural(const MiniStep * pMiniStep) const;

	const Setting * setting(int i) const;

private:
	void preprocessEgo(int ego);
	void preprocessEgo(const Function * pFunction, int ego);
	void calculatePermissibleChanges();
	void calculateTieFlipContributions();
	void calculateTieFlipProbabilities();
	void accumulateScores(int alter) const;
	void accumulateDerivatives() const;
	void copyChangeContributions(MiniStep * pMiniStep) const;
	void checkAlterAgreement(int alter);
	void addAlterAgreementScores(bool accept);
	void accumulateSymmetricModelScores(int alter, bool accept);
// NU	void accumulateRateScores(double tau, const DependentVariable *
//		pSelectedVariable, int selectedActor);
	void calculateSymmetricTieFlipContributions(int alter, int sub);
	void calculateSymmetricTieFlipProbabilities(int alter, int sub, bool aagree);
	bool calculateModelTypeBProbabilities();
	bool diagonalMiniStep(int ego, int alter) const;

	void initializeSetting();
	// The current state of the network
	Network * lpNetwork;

	// The observed data for this network variable
	NetworkLongitudinalData * lpData;

	// The set of actors acting as tie senders
	const SimulationActorSet * lpSenders;

	// The set of actors acting as tie receivers
	const SimulationActorSet * lpReceivers;

	// The number of structural tie variables to active alters per each actor.
	int * lactiveStructuralTieCount {};

	// The current ego in the method makeChange
	int lego {};

	// Indicates if a tie flip to a certain actor is permitted.
	bool * lpermitted {};

	// A two-dimensional array of tie flip contributions to effects, where
	// rows correspond to alters and columns correspond to effects in the
	// evaluation function.

	double ** levaluationEffectContribution {};

	// A two-dimensional array of tie flip contributions to effects, where
	// rows correspond to alters and columns correspond to effects in the
	// endowment function.

	double ** lendowmentEffectContribution {};

	// A two-dimensional array of tie flip contributions to effects, where
	// rows correspond to alters and columns correspond to effects in the
	// tie creation function.

	double ** lcreationEffectContribution {};

	// Selection probability per each alter
	double * lprobabilities {};

	// The cache object for repeated access of various structural properties
	// of this network during the simulation.

	NetworkCache * lpNetworkCache;

	// A vector of filters that may decide that some tie flips from an
	// ego are not permitted.

	std::vector<PermittedChangeFilter *> lpermittedChangeFilters;

	// Vectors of tie flip contributions to effects for models with
	// cooperation between actor and alter.

	double ** lsymmetricEvaluationEffectContribution;
	double ** lsymmetricEndowmentEffectContribution;
	double ** lsymmetricCreationEffectContribution;

	// Probabilities for models with cooperation between actor and alter:
    double lsymmetricProbabilities[2] {};

	// The current alter for two-mode models 
	// with cooperation between actor and alter;
	// and the previous alter for one-mode twostep models:
	int lalter {};
	
	// the probability for taking the double step in option DOUBLESTEP:	
	double lnetworkDoubleStepProb {};
	
	// being one of the option DOUBLESTEP options:	
	bool lnetworkModelTypeDoubleStep {};	
	
	// The probability of selecting the current alter for models with
	// cooperation between actor and alter:
	double lalterProbability {};

	// The probability of making the proposed change for models with
	// cooperation between actor and alter:
	double lsymmetricProbability {};

	ITieIterator* lsettingFlips;

	// whether this is a one mode network or not

	bool loneMode {};

	// the model type:
	NetworkModelType lnetworkModelType;
};


// ----------------------------------------------------------------------------
// Section: Inline accessors
// ----------------------------------------------------------------------------

/**
 * Returns the current ego making the change of the network.
 */
inline int NetworkVariable::ego() const
{
	return this->lego;
}

}

#endif /*NETWORKVARIABLE_H_*/
