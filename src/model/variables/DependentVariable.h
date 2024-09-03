/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DependentVariable.h
 *
 * Description: This file contains the definition of the
 * DependentVariable class.
 *****************************************************************************/

#ifndef DEPENDENTVARIABLE_H_
#define DEPENDENTVARIABLE_H_

#include <map>
#include <vector>
#include <string>
#include "utils/NamedObject.h"
#include "model/Function.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enumerations
// ----------------------------------------------------------------------------

/**
 * This enumeration defines the possible types of model for symmetric,
 * undirected networks;
 * and the possible types of model for behavioral variables.
 */
	enum NetworkModelType { NOTUSED, NORMAL, AFORCE, AAGREE, BFORCE, BAGREE, BJOINT,
		DOUBLESTEP25, DOUBLESTEP50, DOUBLESTEP75, DOUBLESTEP100, NETCONTEMP };
	enum BehaviorModelType { OUTOFUSE, RESTRICT, ABSORB };

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;
class ChangingCovariate;
class NetworkVariable;
class BehaviorVariable;
class EffectValueTable;
class EpochSimulation;
class ActorSet;
class SimulationActorSet;
class LongitudinalData;
class BehaviorLongitudinalData;
class NetworkLongitudinalData;
class Network;
class EffectInfo;
class StructuralRateEffect;
class DiffusionRateEffect;
class MiniStep;
class Setting;


// ----------------------------------------------------------------------------
// Section: DependentVariable class
// ----------------------------------------------------------------------------

/**
 * This class represents a certain dependent variable. It is the base class of
 * NetworkVariable and BehaviorVariable.
 * The class stores the current state of the variable and provides methods
 * supporting simulations of actor-based models.
 */
class DependentVariable : public NamedObject
{
public:
	DependentVariable(std::string name,
		const ActorSet * pActorSet,
		EpochSimulation * pSimulation);
	virtual ~DependentVariable();

	void initializeRateFunction();
	void initializeEvaluationFunction();
	void initializeEndowmentFunction();
	void initializeCreationFunction();

	inline const SimulationActorSet * pActorSet() const;
	int n() const;
	virtual int m() const = 0;
	virtual LongitudinalData * pData() const = 0;
	int id() const;
	virtual bool networkVariable() const;
	virtual bool behaviorVariable() const;
	virtual bool symmetric() const;
	virtual bool constrained() const;
	virtual NetworkModelType networkModelType() const;
	virtual BehaviorModelType behaviorModelType() const;
	virtual bool networkModelTypeB() const;
	virtual bool networkModelTypeDoubleStep() const;
	virtual double networkDoubleStepProb() const;
	virtual int alter() const;

	inline const Function * pEvaluationFunction() const;
	inline const Function * pEndowmentFunction() const;
	inline const Function * pCreationFunction() const;

	virtual void initialize(int period);
	inline int period() const;
	virtual bool canMakeChange(int actor) const;
	virtual void makeChange(int actor) = 0;
	bool successfulChange() const;

	virtual void actOnJoiner(const SimulationActorSet * pActorSet,
		int actor);
	virtual void actOnLeaver(const SimulationActorSet * pActorSet,
		int actor);
	virtual void setLeaverBack(const SimulationActorSet * pActorSet,
		int actor) = 0;

	void calculateRates();
	double totalRate() const;
	double nonSettingsRate() const;
	double rate(int actor) const;
	inline double basicRate() const;
	void updateBasicRate(int period);

	int simulatedDistance() const;

	void accumulateRateScores(double tau,
		const DependentVariable * pSelectedVariable = 0,
		int selectedActor = 0);
	void accumulateRateScores(double tau,
		const DependentVariable * pSelectedVariable,
		int selectedActor, int alter);
	double basicRateScore() const;
	double settingRateScore(std::string setting) const;
	double constantCovariateScore(const ConstantCovariate * pCovariate) const;
	double changingCovariateScore(const ChangingCovariate * pCovariate) const;
	double behaviorVariableScore(const BehaviorVariable * pBehavior) const;
	double outDegreeScore(const NetworkVariable * pNetwork) const;
	double inDegreeScore(const NetworkVariable * pNetwork) const;
	double reciprocalDegreeScore(const NetworkVariable * pNetwork) const;
	double inverseOutDegreeScore(const NetworkVariable * pNetwork) const;
	double logOutDegreeScore(const NetworkVariable * pNetwork) const;
	double inverseInDegreeScore(const NetworkVariable * pNetwork) const;
	double logInDegreeScore(const NetworkVariable * pNetwork) const;
	double inversereciprocalDegreeScore(const NetworkVariable * pNetwork) const;
	double logreciprocalDegreeScore(const NetworkVariable * pNetwork) const;

	// Diffusion effects

	std::map<const EffectInfo *, double> ldiffusionscores;
	std::map<const EffectInfo *, double> ldiffusionsumterms;
	double calculateDiffusionRateEffect(
		const BehaviorVariable * pBehaviorVariable,
		const Network * pNetwork,
		int i, std::string effectName,
		int internalEffectParameter);
	double calculateDiffusionRateEffect(
		const BehaviorVariable * pBehaviorVariable,
		const Network * pNetwork,
		int i, std::string effectName,
		int internalEffectParameter,
		const ConstantCovariate * pConstantCovariate,
		const ChangingCovariate * pChangingCovariate);
// note: calculateDiffusionRateEffect is also a function in StatisticCalculator (almost the same...)

	// Maximum likelihood related

	/**
	 * Calculates the probability of the given ministep assuming that the
	 * ego of the ministep will change this variable.
	 */
	virtual double probability(MiniStep * pMiniStep) = 0;

	virtual bool validMiniStep(const MiniStep * pMiniStep,
		bool checkUpOnlyDownOnlyConditions = true) const;

	void updateEffectParameters();

	/**
	 * Returns if the observed value for the option of the given ministep
	 * is missing at either end of the period.
	 */
	virtual bool missing(const MiniStep * pMiniStep) const = 0;

	/**
	 * Returns if the given ministep is structurally determined for the period
	 */
	virtual bool structural(const MiniStep * pMiniStep) const = 0;

	/**
	 * Generates a random ministep for the given ego.
	 */
	virtual MiniStep * randomMiniStep(int ego) = 0;

	void calculateMaximumLikelihoodRateScores(int activeMiniStepCount);
	void calculateMaximumLikelihoodRateDerivatives(int activeMiniStepCount);
	double basicRateDerivative() const;

	void incrementAcceptances(int stepType);
	void incrementRejections(int stepType);
	void incrementAborts(int stepType);
	int acceptances(int stepType) const;
	int rejections(int stepType) const;
	int aborts(int stepType) const;

	int numberSettings() const;

protected:
	inline EpochSimulation * pSimulation() const;
	void simulatedDistance(int distance);
	void invalidateRates();
	void successfulChange(bool success);
	int stepType() const;
	void getStepType();
	double settingRate() const;
	// A two-dimensional array of tie flip and behavior change contributions to effects.
	std::map<const EffectInfo *, std::vector<double> > * lpChangeContribution;

	Setting** lsettings;
private:
	void initializeFunction(Function * pFunction,
		const std::vector<EffectInfo *> & rEffects) const;

	bool constantRates() const;
	double calculateRate(int i);
	double structuralRate(int i) const;
	double diffusionRate(int i) const;
	double behaviorVariableRate(int i) const;
	void updateCovariateRates();
	void calculateScoreSumTerms();

	// A simulation of the actor-based model, which owns this variable
	EpochSimulation * lpSimulation;

	// The underlying set of actors
	const SimulationActorSet * lpActorSet;

	// The current period (in [0, observations - 2])
	int lperiod {};

	// The total rate of change summed over all actors
	double ltotalRate {};

	// The rate of change factor excepting settingsrate
	double lnonSettingsRate {};

	// The rate of change for each actor
	double * lrate {};

	// The basic rate parameter for the current period
	double lbasicRate {};

	// The setting rate parameters for the current period. Order matches the
	// data object. Only for network variables.
	// double * lsettingRates {};

	// The scaled setting rate parameters for the current period.
	// Order matches the data object. Only for network variables.
	double * lsettingProbs {};

	// The number of settings for this variable.
	// Only non zero for network variables.
	int lnumberSettings {};

	// The type of step in the setting context. -1 if using basic rate or
	// universal setting;
	int lstepType {};

	// The covariate-based component of the rate function per each actor
	double * lcovariateRates {};

	// Parameters for rate effects depending on constant covariates
	std::map<const ConstantCovariate *, double> lconstantCovariateParameters;

	// Parameters for rate effects depending on changing covariates
	std::map<const ChangingCovariate *, double> lchangingCovariateParameters;

	// Parameters for rate effects depending on behavior variables
	std::map<const BehaviorVariable *, double> lbehaviorVariableParameters;

	// The structural rate effects. Currently, there are four types of
	// structural rate effects, namely, the out-degree, in-degree,
	// reciprocal degree, and inverse out-degree effects.

	std::vector<StructuralRateEffect *> lstructuralRateEffects;

	// The diffusion rate effects.

	std::vector<DiffusionRateEffect *> ldiffusionRateEffects;

	// The evaluation function for this variable
	Function * lpEvaluationFunction;

	// The endowment function for this variable
	Function * lpEndowmentFunction;

	// The creation function for this variable
	Function * lpCreationFunction;

	// The distance of this variable from the observed data at the beginning
	// of the current period

	int lsimulatedDistance {};

	// The score for the basic rate parameter for this variable for this period
	double lbasicRateScore {};

	// The derivative for the basic rate parameter for this variable for
	// this period
	double lbasicRateDerivative {};

	// The scores for the setting basic rate parameters for this variable
	// for this period. Only for network variables.
	std::map<std::string, double> lsettingRateScores;

	// Scores for rate effects depending on constant covariates
	std::map<const ConstantCovariate *, double> lconstantCovariateScores;

	// Scores for rate effects depending on changing covariates
	std::map<const ChangingCovariate *, double> lchangingCovariateScores;

	// Scores for rate effects depending on behavior variables
	std::map<const BehaviorVariable *, double> lbehaviorVariableScores;

	// Scores for rate effects depending on out degree
	std::map<const NetworkVariable *, double> loutDegreeScores;

	// Scores for rate effects depending on in degree
	std::map<const NetworkVariable *, double> linDegreeScores;

	// Scores for rate effects depending on reciprocal degree
	std::map<const NetworkVariable *, double> lreciprocalDegreeScores;

	// Scores for rate effects depending on inverse degree
	std::map<const NetworkVariable *, double> linverseOutDegreeScores;

	// Scores for rate effects depending on log degree
	std::map<const NetworkVariable *, double> llogOutDegreeScores;

	// Scores for rate effects depending on inverse in degree
	std::map<const NetworkVariable *, double> linverseInDegreeScores;

	// Scores for rate effects depending on log in degree
	std::map<const NetworkVariable *, double> llogInDegreeScores;

	// Scores for rate effects depending on inverse reciprocal degree
	std::map<const NetworkVariable *, double> linversereciprocalDegreeScores;

	// Scores for rate effects depending on log reciprocal degree
	std::map<const NetworkVariable *, double> llogreciprocalDegreeScores;

	// Sum term for scores for rate effects depending on constant covariates
	std::map<const ConstantCovariate *, double> lconstantCovariateSumTerm;

	// Sum term for scores for rate effects depending on changing covariates
	std::map<const ChangingCovariate *, double> lchangingCovariateSumTerm;

	// Sum term for scores for rate effects depending on behavior variables
	std::map<const BehaviorVariable *, double> lbehaviorVariableSumTerm;

	// Sum term for scores for rate effects depending on out degree
	std::map<const NetworkVariable *, double> loutDegreeSumTerm;

	// Sum term for scores for rate effects depending on in degree
	std::map<const NetworkVariable *, double> linDegreeSumTerm;

	// Sum term for scores for rate effects depending on reciprocal degree
	std::map<const NetworkVariable *, double> lreciprocalDegreeSumTerm;

	// Sum term for scores for rate effects depending on inverse degree
	std::map<const NetworkVariable *, double> linverseOutDegreeSumTerm;

	// Sum term for scores for rate effects depending on log degree
	std::map<const NetworkVariable *, double> llogOutDegreeSumTerm;

	// Sum term for scores for rate effects depending on inverse degree
	std::map<const NetworkVariable *, double> linverseInDegreeSumTerm;

	// Sum term for scores for rate effects depending on log degree
	std::map<const NetworkVariable *, double> llogInDegreeSumTerm;

	// Sum term for scores for rate effects depending on inverse reciprocal degree
	std::map<const NetworkVariable *, double> linversereciprocalDegreeSumTerm;

	// Sum term for scores for rate effects depending on log reciprocal degree
	std::map<const NetworkVariable *, double> llogreciprocalDegreeSumTerm;

	// Sum term for model B scores for rate effects depending on constant
	// covariates
	std::map<const ConstantCovariate *, double> lconstantCovariateModelBSumTerm;

	// Sum term for model B scores for rate effects depending on changing
	// covariates
	std::map<const ChangingCovariate *, double> lchangingCovariateModelBSumTerm;

	// Sum term for model B scores for rate effects depending on behavior
	// variables
	std::map<const BehaviorVariable *, double> lbehaviorVariableModelBSumTerm;

	// Sum term for model B scores for rate effects depending on out degree
	std::map<const NetworkVariable *, double> loutDegreeModelBSumTerm;

	// Sum term for model B scores for rate effects depending on inverse degree
	std::map<const NetworkVariable *, double> linverseOutDegreeModelBSumTerm;

	// Sum term for model B scores for rate effects depending on log degree
	std::map<const NetworkVariable *, double> llogOutDegreeModelBSumTerm;

	// Sum term for model B scores for rate effects depending on inverse in degree
	std::map<const NetworkVariable *, double> linverseInDegreeModelBSumTerm;

	// Sum term for model B scores for rate effects depending on log in degree
	std::map<const NetworkVariable *, double> llogInDegreeModelBSumTerm;

	// Indicates if the rates are valid and shouldn't be recalculated
	// provided that the rates are constant during the period.

	int lvalidRates {};

	// flag to indicate we gave up on a step due to uponly and other filters
	bool lsuccessfulChange {};

	std::vector <int> lacceptances {};
	std::vector <int> lrejections {};
	std::vector <int> laborts {};

};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the actor-based model owning this variable.
 */
EpochSimulation * DependentVariable::pSimulation() const
{
	return this->lpSimulation;
}


/**
 * Returns the set of actors underlying this dependent variable.
 */
const SimulationActorSet * DependentVariable::pActorSet() const
{
	return this->lpActorSet;
}


/**
 * Returns the evaluation function of this variable.
 */
const Function * DependentVariable::pEvaluationFunction() const
{
	return this->lpEvaluationFunction;
}


/**
 * Returns the endowment function of this variable.
 */
const Function * DependentVariable::pEndowmentFunction() const
{
	return this->lpEndowmentFunction;
}


/**
 * Returns the tie creation function of this variable.
 */
const Function * DependentVariable::pCreationFunction() const
{
	return this->lpCreationFunction;
}


/**
 * Returns the index of the current period.
 */
int DependentVariable::period() const
{
	return this->lperiod;
}


/**
 * Returns the basic rate parameter for the current period.
 */
double DependentVariable::basicRate() const
{
	return this->lbasicRate;
}

}

#endif /*DEPENDENTVARIABLE_H_*/
