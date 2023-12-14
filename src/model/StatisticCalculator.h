/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: StatisticCalculator.h
 *
 * Description: This file contains the definition of the
 * StatisticCalculator class.
 *****************************************************************************/

#ifndef STATISTICCALCULATOR_H_
#define STATISTICCALCULATOR_H_

#include <map>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Data;
class Model;
class State;
class EffectInfo;
class LongitudinalData;
class NetworkLongitudinalData;
class BehaviorLongitudinalData;
class ContinuousLongitudinalData;
class Network;
class EffectFactory;
class ConstantCovariate;
class ChangingCovariate;

// ----------------------------------------------------------------------------
// Section: StatisticCalculator class
// ----------------------------------------------------------------------------

/**
 * Provides means for calculating statistics corresponding to effects of a
 * certain model. In order to calculate the statistics for effects of a
 * model, one has to provide the underlying observed data, the period under
 * consideration, and the current state of all dependent variables. This
 * information should be given at the construction of the statistic calculator,
 * after which individual statistics can be queried. Example:
 *
 * StatisticCalculator calculator(pData, pModel, pState, period);
 * double statistic1 = calculator.statistic(pEffectInfo1);
 * double statistic2 = calculator.statistic(pEffectInfo2);
 */
class StatisticCalculator
{
public:
	StatisticCalculator(const Data * pData,
		const Model * pModel,
		State * pState,
		int period);
	StatisticCalculator(const Data * pData,
		const Model * pModel, State * pState,
		int period, bool returnActorStatistics);
	StatisticCalculator(const Data * pData,
		const Model * pModel, State * pState,
		int period, bool returnActorStatistics, bool returnStaticChangeContributions);
	virtual ~StatisticCalculator();

	double statistic(EffectInfo * pEffectInfo) const;
	std::vector<double *> staticChangeContributions(EffectInfo * pEffect) const;
	double * actorStatistics(EffectInfo * pEffect) const;
	int distance(LongitudinalData * pData, int period) const;
	double distance(ContinuousLongitudinalData * pData, int period) const;
	double totalDistance(int period) const;
	int settingDistance(LongitudinalData * pData, std::string setting,
		int period) const;

private:
	void calculateStatisticsInitNetwork(NetworkLongitudinalData * pNetwork);
	void calculateStatistics();
	void calculateNetworkRateStatistics(
		NetworkLongitudinalData * pNetworkData);
	void calculateNetworkGMMStatistics(NetworkLongitudinalData * pNetworkData);
	void calculateNetworkEvaluationStatistics(
		NetworkLongitudinalData * pNetworkData);
	void calculateNetworkEndowmentStatistics(
		NetworkLongitudinalData * pNetworkData);
	void calculateNetworkCreationStatistics(
		NetworkLongitudinalData * pNetworkData);
	void calculateBehaviorStatistics(BehaviorLongitudinalData * pBehaviorData);
	void calculateBehaviorGMMStatistics(
			BehaviorLongitudinalData * pBehaviorData);
	void calculateBehaviorRateStatistics(BehaviorLongitudinalData * pBehaviorData);
	void calculateContinuousStatistics(ContinuousLongitudinalData * pContinuousData);
	void calculateContinuousRateStatistics(ContinuousLongitudinalData * pContinuousData);
	// Functions to calculate value of diffusion rate effect
	// note: calculateDiffusionRateEffect is also a function in DependentVariable (almost the same...)
	double calculateDiffusionRateEffect(BehaviorLongitudinalData *
		pBehaviorData, const Network * pStructural, int i,
		std::string effectName, int internalEffectParameter);
	double calculateDiffusionRateEffect(BehaviorLongitudinalData *
		pBehaviorData, const Network * pStructural,
		const ConstantCovariate * pConstantCovariate,
		const ChangingCovariate * pChangingCovariate,
		int i, std::string effectName, int internalEffectParameter);

	// The data to be used for calculating the statistics
	const Data * lpData;

	// The model containing the effects whose statistics are to be calculated
	const Model * lpModel;

	// The current state of the dependent variables
	State * lpState;

	// The period of the evolution we are interested in
	int lperiod {};

	// indicates whether actor statistics are needed
	bool lneedActorStatistics {};
	
	// indicates whether change contributions are needed
	bool lcountStaticChangeContributions {};

	// The resulting map of statistic values
	std::map<EffectInfo *, double> lstatistics;

	// The resulting map of actor statistic values
	std::map<EffectInfo *, double * > lactorStatistics;

	// The change contributions of all effects
	std::map<EffectInfo *, std::vector<double *> > lstaticChangeContributions;

	// Array of simulated distances per variable
	std::map<LongitudinalData *, int *> ldistances;

	// Array of simulated distances per continuous variable
	std::map<ContinuousLongitudinalData *, double *> lcontinuousDistances;
	
	// Array of simulated distances per setting per network variable
	std::map<LongitudinalData *, std::map<std::string, int *> > lsettingDistances;

	State * lpPredictorState;

	State * lpStateLessMissingsEtc;

	void calcDifferences(
			NetworkLongitudinalData * const pNetworkData,
			const Network* const pDifference);
};

}

#endif /*STATISTICCALCULATOR_H_*/
