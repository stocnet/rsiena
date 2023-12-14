/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Model.h
 *
 * Description: This file contains the definition of the Model class.
 *****************************************************************************/

#ifndef MODEL_H_
#define MODEL_H_

#include <map>
#include <vector>
#include <string>

namespace siena
{
// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Data;
class LongitudinalData;
class NetworkLongitudinalData;
class Effect;
class StructuralRateEffect;
class Function;
class EffectInfo;
class Chain;


// ----------------------------------------------------------------------------
// Section: Model class
// ----------------------------------------------------------------------------

/**
 * This class defines the actor-based models for longitudinal network and
 * behavioral data.
 */
class Model
{
public:
	Model();
	virtual ~Model();

	// Basic rate and scale effects

	void basicRateParameter(LongitudinalData * pDependentVariableData,
		int period,
		double value);
	double basicRateParameter(LongitudinalData * pDependentVariableData,
		int period) const;
	void basicScaleParameter(int period, double value);
	double basicScaleParameter(int period) const;		

	// Setting rate effects

	void settingRateParameter(NetworkLongitudinalData * pNetworkData,
		std::string setting,
		int period,
		double value);
	double settingRateParameter(NetworkLongitudinalData * pNetworkData,
		std::string setting,
		int period) const;
	int numberOfSettings(NetworkLongitudinalData * pNetworkData) const;

	// Other effects

	EffectInfo * addEffect(std::string variableName,
		std::string effectName,
		std::string effectType,
		double parameter,
		double internalEffectParameter = 0,
		std::string interactionName1 = "",
		std::string interactionName2 = "",
		std::string rateType = "");
	EffectInfo * addInteractionEffect(std::string variableName,
		std::string effectName,
		std::string effectType,
		double parameter,
		const EffectInfo * pEffect1,
		const EffectInfo * pEffect2,
		const EffectInfo * pEffect3 = 0);

	bool gmmModel() const;

	const std::vector<EffectInfo *> & rRateEffects(std::string variableName) const;
	const std::vector<EffectInfo *> & rGMMEffects(std::string variableName) const;
	const std::vector<EffectInfo *> & rEvaluationEffects(std::string variableName) const;
	const std::vector<EffectInfo *> & rEndowmentEffects(std::string variableName) const;
	const std::vector<EffectInfo *> & rCreationEffects(std::string variableName) const;

	void chainStore(const Chain& chain, int periodFromStart);
	std::vector <Chain *> & rChainStore( int periodFromStart);
	void clearChainStore(int keep, int periodFromStart);
	void setupChainStore(int numberOfPeriods);
	void deleteLastChainStore(int periodFromStart);
	void numberOfPeriods(int numberOfPeriods);
	int numberOfPeriods();

	// Various flags

	void conditional(bool flag);
	bool conditional() const;

	int targetChange(const Data * pData, int period) const;
	void targetChange(const Data * pData, int period, int change);

	std::string conditionalDependentVariable() const;
	void conditionalDependentVariable(std::string variableName);

	void needChain(bool flag);
	bool needChain() const;

	void needScores(bool flag);
	bool needScores() const;

	void needDerivatives(bool flag);
	bool needDerivatives() const;

	void parallelRun(bool flag);
	bool parallelRun() const;

	void needChangeContributions(bool flag);
	bool needChangeContributions() const;

	// various stores for ML

	void numberMLSteps(int value);
	int numberMLSteps() const;

	void maximumPermutationLength(double value);
	double maximumPermutationLength() const;

	void minimumPermutationLength(double value);
	double minimumPermutationLength() const;

	void initialPermutationLength(double value);
	double initialPermutationLength() const;

	void initializeCurrentPermutationLength();
	double currentPermutationLength(int period) const;
	void currentPermutationLength(int period, double value);

	void insertDiagonalProbability(double probability);
	double insertDiagonalProbability() const;

	void cancelDiagonalProbability(double probability);
	double cancelDiagonalProbability() const;

	void permuteProbability(double probability);
	double permuteProbability() const;

	void insertPermuteProbability(double probability);
	double insertPermuteProbability() const;

	void deletePermuteProbability(double probability);
	double deletePermuteProbability() const;

	void insertRandomMissingProbability(double probability);
	double insertRandomMissingProbability() const;

	void deleteRandomMissingProbability(double probability);
	double deleteRandomMissingProbability() const;

	void missingNetworkProbability(double probability);
	double missingNetworkProbability(int periodFromStart) const;

	void missingBehaviorProbability(double probability);
	double missingBehaviorProbability(int periodFromStart) const;

	void normalizeSettingRates(bool normalize);
	bool normalizeSettingRates() const;

	// localML
	void localML(bool flag);
	bool localML() const;

	// simple rates flag for ML
	bool simpleRates() const;
	void simpleRates(bool simpleRates);

private:
	// Indicates if conditional simulation has to be carried out
	bool lconditional {};

	//! True if any dependent variable has effects in the gmm objective
	//! function.
	bool lGMMModel {};

	// name of conditional dependent variable
	std::string lconditionalDependentVariable {};

	// Targets for conditional dependent variable per each data object
	// and period.

	std::map<const Data *, int *> ltargetChanges;

	// An array of doubles per each longitudinal data object storing
	// the basic rate parameters for all periods

	std::map<const LongitudinalData *, double *> lbasicRateParameters;

	// An array of doubles for some network longitudinal data objects storing
	// the basic rate parameters for all periods by setting.

	std::map<const NetworkLongitudinalData *, std::map<std::string, double *> >
		lsettingRateParameters;

	// An array of doubles storing the scale parameters for all periods
	// (part of the SDE model for continuous variables)
	double * lbasicScaleParameters {};

	// A vector of effects other than the basic rate effects.
	std::vector<EffectInfo *> leffects;

	// A vector of rate effects (except the basic rate effects) per variable
	std::map<std::string, std::vector<EffectInfo *> > lrateEffects;

	// A vector of pointers to GMM effects per variable
	std::map<std::string, std::vector<EffectInfo *> > lgmmEffects;

	// A vector of pointers to evaluation effects per variable
	std::map<std::string, std::vector<EffectInfo *> > levaluationEffects;

	// A vector of pointers to endowment effects per variable
	std::map<std::string, std::vector<EffectInfo *> > lendowmentEffects;

	// A vector of pointers to creation effects per variable
	std::map<std::string, std::vector<EffectInfo *> > lcreationEffects;

	// A dummy vector of effect infos in case we need a reference to
	// non-existent vectors

	const std::vector<EffectInfo *> lemptyEffectVector;

	// indicates whether we need to keep a chain of ministeps
	bool lneedChain {};

	// indicates whether we need to accumulate scores in this iteration
	bool lneedScores {};

	// indicates whether we need to accumulate derivatives for ML in
	// this iteration
	bool lneedDerivatives {};

	// indicates whether we need to store the change contributions on
	// the ministeps
	bool lneedChangeContributions {};

	// indicates whether we need to match Siena3 in use of random variables
	// and score calculations
	bool lparallelRun {};

	//indicates whether change contributions are needed
	bool lneedChangeContributions2 {};

	// number of steps in a run for ML
	int lnumberMLSteps {};

	// maximum length of permuted interval
	double lmaximumPermutationLength {};

	// minimum length of permuted interval
	double lminimumPermutationLength {};

	// initial length of permuted interval
	double linitialPermutationLength {};

	// current length of permuted interval: varies by period
	std::vector <double> lcurrentPermutationLength;

	// probabilities of the different ML steps
	double linsertDiagonalProbability {};
	double lcancelDiagonalProbability {};
	double lpermuteProbability {};
	double linsertPermuteProbability {};
	double ldeletePermuteProbability {};
	double linsertRandomMissingProbability {};
	double ldeleteRandomMissingProbability {};

	// localML
	bool llocalML {};

	bool lsimpleRates {};

	std::vector <double> lmissingNetworkProbability;
	std::vector <double> lmissingBehaviorProbability;

	// chain storage: vector of chains for each period for each set of samples
	// lchainStore[i] is set of entries for periodFromStart <i>,
	// which incorporates both the group and period.
	std::vector <std::vector <Chain *> > lchainStore;

	int lnumberOfPeriods {};
	bool lnormalizeSettingsRates {};
};

}

#endif /*MODEL_H_*/
