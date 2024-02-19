/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: StatisticCalculator.cpp
 *
 * Description: This file contains the implementation of the
 * StatisticCalculator class.
 *****************************************************************************/

#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <R_ext/Arith.h>
#include <R_ext/Print.h>

#include "StatisticCalculator.h"
#include "data/Data.h"
#include "data/ActorSet.h"
#include "network/NetworkUtils.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ContinuousLongitudinalData.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingDyadicCovariate.h"
#include "model/Model.h"
#include "model/State.h"
#include "model/Function.h"
#include "model/EffectInfo.h"
#include "model/effects/EffectFactory.h"
#include "model/effects/NetworkEffect.h"
#include "model/effects/BehaviorEffect.h"
#include "model/effects/ContinuousEffect.h"
#include "model/EpochSimulation.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/tables/Cache.h"
#include "network/IncidentTieIterator.h"
#include "network/layers/DistanceTwoLayer.h"
#include "network/iterators/UnionTieIterator.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction, destruction
// ----------------------------------------------------------------------------

/**
 * Constructor.
 * @param[in] pData the observed data
 * @param[in] pModel the model whose effect statistics are to be calculated
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period under consideration
 */
StatisticCalculator::StatisticCalculator(const Data * pData,
	const Model * pModel, State * pState, int period)
{
	this->lpData = pData;
	this->lpModel = pModel;
	this->lpState = pState;
	this->lperiod = period;
	this->lpPredictorState = new State();
	this->lpStateLessMissingsEtc = new State();
	this->lneedActorStatistics = 0;
	this->lcountStaticChangeContributions = 0;

	this->calculateStatistics();
}

/**
 * Constructor.
 * @param[in] pData the observed data
 * @param[in] pModel the model whose effect statistics are to be calculated
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period under consideration
 * @param[in] returnActorStatistics whether individual actor statistics should be returned
 */
StatisticCalculator::StatisticCalculator(const Data * pData,
		const Model * pModel, State * pState, int period,
		bool returnActorStatistics)
{
	this->lpData = pData;
	this->lpModel = pModel;
	this->lpState = pState;
	this->lperiod = period;
	this->lpPredictorState = new State();
	this->lpStateLessMissingsEtc = new State();
	this->lneedActorStatistics = returnActorStatistics;
	this->lcountStaticChangeContributions = 0;

	this->calculateStatistics();
}

/**
 * Constructor.
 * @param[in] pData the observed data
 * @param[in] pModel the model whose effect statistics are to be calculated
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period under consideration
 * @param[in] returnActorStatistics whether individual actor statistics should be returned
 * @param[in] returnStaticChangeContributions whether contributions of effects on possible next tie flips or behavior changes are needed
 */
StatisticCalculator::StatisticCalculator(const Data * pData,
		const Model * pModel, State * pState, int period,
		bool returnActorStatistics, bool returnStaticChangeContributions)
{
	this->lpData = pData;
	this->lpModel = pModel;
	this->lpState = pState;
	this->lperiod = period;
	this->lpPredictorState = new State();
	this->lpStateLessMissingsEtc = new State();
	this->lneedActorStatistics = returnActorStatistics;
	this->lcountStaticChangeContributions = returnStaticChangeContributions;

	this->calculateStatistics();
}

template<typename T, typename A> // type and allocator
static void clear_vector_of_array_pointers(vector<T*, A>& v)
{
	for (typename vector<T*, A>::iterator it = v.begin(); it != v.end(); ++it)
		delete[] (*it);
	v.clear();
}

template <typename K, typename T>
static void clear_map_value_array_pointers(map<K, T*>& m)
{
	for (typename map<K, T*>::iterator it = m.begin(); it != m.end(); ++it)
		delete[] it->second;
	m.clear();
}

template <typename K, typename T>
static void for_each_map_value(map<K, T>& m, void (*fn)(T&))
{
	for (typename map<K, T>::iterator it = m.begin(); it != m.end(); ++it)
		(*fn)(it->second);
}

/**
 * Deallocates this model.
 */
StatisticCalculator::~StatisticCalculator()
{
	clear_map_value_array_pointers(this->ldistances);
	clear_map_value_array_pointers(this->lcontinuousDistances);

	for_each_map_value(this->lsettingDistances, &clear_map_value_array_pointers);
	this->lsettingDistances.clear();

	// this->lpPredictorState->deleteValues(); // now called by the State dtor
	delete this->lpPredictorState;
	this->lpPredictorState = 0;

	// this->lpStateLessMissingsEtc->deleteValues(); // now called by the State dtor
	delete this->lpStateLessMissingsEtc;
	this->lpStateLessMissingsEtc = 0;

	for_each_map_value(this->lstaticChangeContributions,
			&clear_vector_of_array_pointers);
	this->lstaticChangeContributions.clear();

	clear_map_value_array_pointers(this->lactorStatistics);
}

// ----------------------------------------------------------------------------
// Section: Accessors for statistics
// ----------------------------------------------------------------------------

/**
 * Returns the statistic for the given effect.
 */
double StatisticCalculator::statistic(EffectInfo * pEffect) const
{
	map<EffectInfo *, double>::const_iterator iter =
		this->lstatistics.find(pEffect);

	if (iter == this->lstatistics.end())
	{
		throw invalid_argument(
			"Unknown effect: The given effect is not part of the model.");
	}

	return iter->second;
}

/**
 * Returns the tie flip contributions or the behavior change contributions of the given effect.
 */
vector<double *> StatisticCalculator::staticChangeContributions(EffectInfo * pEffect) const
{
	map<EffectInfo *, vector<double *> >::const_iterator iter =
		this->lstaticChangeContributions.find(pEffect);
	if (iter == this->lstaticChangeContributions.end())
	{
		throw invalid_argument(
			"Unknown effect: The given effect is not part of the model.");
	}
	return iter->second;
}

/**
 * Returns the actor statistics of the given effect.
 */
double * StatisticCalculator::actorStatistics(EffectInfo * pEffect) const
{
	map<EffectInfo *, double * >::const_iterator iter =
		this->lactorStatistics.find(pEffect);

	if (iter == this->lactorStatistics.end())
	{
		throw invalid_argument(
			"Unknown effect: The given effect is not part of the model.");
	}
	return iter->second;
}

/**
 * Returns the simulated distance for the given network and period.
 */
int StatisticCalculator::distance(LongitudinalData * pData, int period)
	const
{
	map<LongitudinalData *, int *>::const_iterator iter =
		this->ldistances.find(pData);

	if (iter == this->ldistances.end())
	{
		throw invalid_argument(
			"Unknown effect: The given basic rate is not part of the model.");
	}

	return iter->second[period];
}


/**
 * Returns the simulated distance for the given continuous behavior variable
 * and period.
 */
double StatisticCalculator::distance(ContinuousLongitudinalData * pData, int period)
	const
{
	map<ContinuousLongitudinalData *, double *>::const_iterator iter =
		this->lcontinuousDistances.find(pData);

	if (iter == this->lcontinuousDistances.end())
	{
		throw invalid_argument(
			"Unknown effect: The given scale parameter is not part of the model.");
	}

	return iter->second[period];
}


/**
 * Returns the total simulated distance of all continuous behavior variables
 * for the given period.
 */
double StatisticCalculator::totalDistance(int period) const
{
	double total = 0;

	for (map<ContinuousLongitudinalData *, double *>::const_iterator iter =
			this->lcontinuousDistances.begin();
		 iter != this->lcontinuousDistances.end();
		 iter++)
	{
		total += iter->second[period];
	}

	return total;
}


/**
 * Returns the simulated setting distance for the given network and period.
 */
int StatisticCalculator::settingDistance(LongitudinalData * pData,
	string setting, int period) const
{
	int value = 0;
	map<LongitudinalData *, map<string, int *> >::const_iterator iter=
 		this->lsettingDistances.find(pData);

	if (iter != this->lsettingDistances.end())
	{
		map<string, int *>::const_iterator iter1 =
			iter->second.find(setting);
		value = iter1->second[period];
	}
	else
	{
		throw invalid_argument("Unknown setting: " + setting);
	}

	return value;
}

void StatisticCalculator::calculateStatisticsInitNetwork(NetworkLongitudinalData * pNetworkData)
{
	const Network * pPredictor = pNetworkData->pNetworkLessMissing(this->lperiod);
	this->lpPredictorState->pNetwork(pNetworkData->name(), pPredictor);

	// Duplicate the current network and remove those ties that are missing at
	// either end of the period.
			Network * pNetwork = this->lpState->pNetwork(pNetworkData->name())->clone();

			subtractNetwork(pNetwork, pNetworkData->pMissingTieNetwork(this->lperiod));
			subtractNetwork(pNetwork, pNetworkData->pMissingTieNetwork(this->lperiod + 1));

			// for not-targets, overwrite the current network for values
			// structurally fixed for the next period. (no effect for targets)

			replaceNetwork(pNetwork,
				pNetworkData->pNetwork(this->lperiod + 1),
				pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

			// for targets look backwards and mimic the simulation by carrying
			// forward structural values.

			replaceNetwork(pNetwork,
				pNetworkData->pNetwork(this->lperiod),
				pNetworkData->pStructuralTieNetwork(this->lperiod));

			// NOTE: pass delete responsibility to state
	this->lpStateLessMissingsEtc->pNetwork(pNetworkData->name(), pNetwork);
			// delete pNetwork;
}

/**
 * Calculates the statistics for all effects of the given model. Note that
 * this->lperiod relates to the current period when simulating, but the
 * previous when calculating targets.
 */
void StatisticCalculator::calculateStatistics()
{
	const vector<LongitudinalData *> & rVariables = this->lpData->rDependentVariableData();

	// set up the predictor and currentLessMissingsEtc states of these variables
	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData = dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
		BehaviorLongitudinalData * pBehaviorData = dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);
		ContinuousLongitudinalData * pContinuousData = dynamic_cast<ContinuousLongitudinalData *>(rVariables[i]);

		if (pNetworkData)
		{
			calculateStatisticsInitNetwork(pNetworkData);
		}
		else if (pBehaviorData)
		{
			// create a copy of the start of the period and zero any values missing
			// at (either end?) start of period
			const int * values = pBehaviorData->valuesLessMissingStarts(this->lperiod);
			this->lpPredictorState->behaviorValues(pBehaviorData->name(), values);
		}
		else if (pContinuousData)
		{
			const double * values = pContinuousData->valuesLessMissingStarts(this->lperiod);
			this->lpPredictorState->continuousValues(pContinuousData->name(), values);
		}
		else
		{
			throw domain_error("Unexpected class of dependent variable");
		}
	}

	const vector<LongitudinalData *> & rSimVariables = this->lpData->rSimVariableData();
	for (unsigned i = 0; i < rSimVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData = dynamic_cast<NetworkLongitudinalData *>(rSimVariables[i]);

		if (pNetworkData)
		{
			calculateStatisticsInitNetwork(pNetworkData);
//			Network * pNetwork = this->lpState->pNetwork(pNetworkData->name())->clone();
//			this->lpPredictorState->pNetwork(pNetworkData->name(), pNetwork);
		}
		else
		{
			throw domain_error("Unexpected class of simulated variable");
		}
	}

	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);
		ContinuousLongitudinalData * pContinuousData =
			dynamic_cast<ContinuousLongitudinalData *>(rVariables[i]);

		if (pNetworkData)
		{
			this->calculateNetworkRateStatistics(pNetworkData);
			this->calculateNetworkEvaluationStatistics(pNetworkData);
			this->calculateNetworkEndowmentStatistics(pNetworkData);
			this->calculateNetworkCreationStatistics(pNetworkData);
			this->calculateNetworkGMMStatistics(pNetworkData);
		}
		else if (pBehaviorData)
		{
			this->calculateBehaviorRateStatistics(pBehaviorData);
			this->calculateBehaviorStatistics(pBehaviorData);
			this->calculateBehaviorGMMStatistics(pBehaviorData);
		}
		else if (pContinuousData)
		{
			this->calculateContinuousRateStatistics(pContinuousData);
			this->calculateContinuousStatistics(pContinuousData);
		}
		else
		{
			throw domain_error("Unexpected class of dependent variable");
		}
	}

	// clean up created data not owned by states
	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
		string name = rVariables[i]->name();

		if (pNetworkData)
		{
			const Network * pNetwork = this->lpStateLessMissingsEtc->pNetwork(name);
			delete pNetwork;
		}
	}
}

/**
 * Calculates the statistics for the GMM effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkGMMStatistics(
		NetworkLongitudinalData * pNetworkData) {
	// We want to pass all networks to the effects in a single state,
	// hence we overwrite the network in the predictor state.
	// We do not use the predictor network with effects this network owns.

	string name = pNetworkData->name();
	const Network * pPredictorNetwork = this->lpPredictorState->pNetwork(name);
	const Network * pCurrentLessMissingsEtc =
			this->lpStateLessMissingsEtc->pNetwork(name);
	this->lpPredictorState->pNetwork(name, pCurrentLessMissingsEtc);

	// Loop through the evaluation effects, calculate the statistics, and store
	// them.
	const vector<EffectInfo *> & rEffects = this->lpModel->rGMMEffects(
			pNetworkData->name());
	EffectFactory factory(this->lpData);
	Cache cache;

	for (unsigned i = 0; i < rEffects.size(); i++) {
		EffectInfo * pInfo = rEffects[i];
		NetworkEffect * pEffect = (NetworkEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.
//	pEffect->initialize(this->lpData, this->lpPredictorState,
//		this->lperiod, &cache);
		pEffect->initialize(this->lpData, this->lpPredictorState, this->lpState,
			this->lperiod, &cache);

		this->lstatistics[pInfo] = pEffect->evaluationStatistic();
		delete pEffect;
	}

	// Restore the predictor network
	this->lpPredictorState->pNetwork(name, pPredictorNetwork);
}

void StatisticCalculator::calculateBehaviorGMMStatistics(
		BehaviorLongitudinalData * pBehaviorData) {
	// create a copy of the current state and zero any values missing
	// at either end of period
	const int* currentState = this->lpState->behaviorValues(
			pBehaviorData->name());
	double* currentValues = new double[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++) {
		currentValues[i] = currentState[i] - pBehaviorData->overallMean();
		if (pBehaviorData->missing(this->lperiod, i)
				|| pBehaviorData->missing(this->lperiod + 1, i)) {
			currentValues[i] = 0;
		}
	}
	// Loop through the evaluation effects, calculate the statistics,
	// and store them.
	const vector<EffectInfo *> & rEffects = this->lpModel->rGMMEffects(
			pBehaviorData->name());
	EffectFactory factory(this->lpData);
	Cache cache;
	for (unsigned i = 0; i < rEffects.size(); i++) {
		EffectInfo * pInfo = rEffects[i];
		BehaviorEffect * pEffect = (BehaviorEffect *) factory.createEffect(
				pInfo);
		// Initialize the effect to work with our data and state of variables.
//	pEffect->initialize(this->lpData, this->lpPredictorState,
//			this->lperiod, &cache);
		pEffect->initialize(this->lpData, this->lpPredictorState, this->lpState,
				this->lperiod, &cache);
		if (this->lneedActorStatistics) {
			pair<double, double *> p = pEffect->evaluationStatistic(
					currentValues, this->lneedActorStatistics);
			this->lstatistics[pInfo] = p.first;
			this->lactorStatistics[pInfo] = p.second;
		} else {
			this->lstatistics[pInfo] = pEffect->evaluationStatistic(
					currentValues);
		}
		if (this->lcountStaticChangeContributions)
		{
			int choices = 3;
			vector<double *> egosMap(pBehaviorData->n());
			this->lstaticChangeContributions.insert(make_pair(pInfo, egosMap));
			for (int e = 0; e < pBehaviorData->n(); e++) {
				cache.initialize(e);
				//			pEffect->initialize(this->lpData, this->lpPredictorState,
				//					this->lperiod, &cache);
				pEffect->initialize(this->lpData, this->lpPredictorState, this->lpState,
						this->lperiod, &cache);
				double * contributions = new double[choices];
				this->lstaticChangeContributions[pInfo].at(e) = contributions;
				pEffect->preprocessEgo(e);
				// no change gives no contribution
				this->lstaticChangeContributions[pInfo].at(e)[1] = 0;
				// calculate the contribution of downward change
				if ((currentState[e] > pBehaviorData->min())  // was currentValues
						&& (!pBehaviorData->upOnly(this->lperiod)))
				{
					this->lstaticChangeContributions[pInfo].at(e)[0] =
						pEffect->calculateChangeContribution(e, -1);
				} else {
					this->lstaticChangeContributions[pInfo].at(e)[0] = R_NaN;
				}
				// calculate the contribution of upward change
				if ((currentState[e] < pBehaviorData->max())  // was currentValues
						&& (!pBehaviorData->downOnly(this->lperiod)))
				{
					this->lstaticChangeContributions[pInfo].at(e)[2] =
						pEffect->calculateChangeContribution(e, 1);
				} else {
					this->lstaticChangeContributions[pInfo].at(e)[2] = R_NaN;
				}
			}
		}
		delete pEffect;
	}
	delete[] currentValues;
}

/**
 * Calculates the statistics for the evaluation effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkEvaluationStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	// We want to pass all networks to the effects in a single state,
	// hence we overwrite the network in the predictor state.
	// We do not use the predictor network with effects this network owns.

	string name = pNetworkData->name();
	const Network * pPredictorNetwork = this->lpPredictorState->pNetwork(name);
	const Network * pCurrentLessMissingsEtc =
		this->lpStateLessMissingsEtc->pNetwork(name);
	this->lpPredictorState->pNetwork(name, pCurrentLessMissingsEtc);

	// Loop through the evaluation effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEffects =
		this->lpModel->rEvaluationEffects(pNetworkData->name());

	EffectFactory factory(this->lpData);
	Cache cache;

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		NetworkEffect * pEffect = (NetworkEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData, this->lpPredictorState,
			this->lperiod, &cache);

		if(this->lneedActorStatistics)
		{
			pair<double, double * > p = pEffect->evaluationStatistic(this->lneedActorStatistics);
			this->lstatistics[pInfo] = p.first;
			this->lactorStatistics[pInfo] = p.second;
		}
		else
		{
			this->lstatistics[pInfo] = pEffect->evaluationStatistic();
		}

// I think the following is used only for sienaRI; note the "static", this is used for observations only.
		if (this->lcountStaticChangeContributions)
		{			
//	Rprintf(" this->lcountStaticChangeContributions in calculateNetworkEvaluationStatistics \n");
			int egos = pCurrentLessMissingsEtc->n();
			int alters = pCurrentLessMissingsEtc->m();
			vector<double *> egosMap(egos);
			this->lstaticChangeContributions.insert(make_pair(pInfo,egosMap));
			for (int e = 0; e < egos ; e++)
			{
				cache.initialize(e);
				pEffect->initialize(this->lpData, this->lpPredictorState,
						this->lperiod, &cache);
				double * contributions = new double[egos];
					this->lstaticChangeContributions[pInfo].at(e) = contributions;
					pEffect->preprocessEgo(e);
					// TODO determine permissible changes
					// (see: NetworkVariable::calculatePermissibleChanges())
					// TODO allow also bipartite; this requires knowledge of number of alters.
//earlier:					for (int a = 0; a < alters ; a++)
					for (int a = 0; a < egos ; a++)
					{
						if ((a == e) && (pNetworkData->oneModeNetwork()))
						{
							this->lstaticChangeContributions[pInfo].at(e)[a] = 0;
						}
						else if (a == alters)
						{
							this->lstaticChangeContributions[pInfo].at(e)[a] = 0;
						}
						else if (a > alters)
						{
							this->lstaticChangeContributions[pInfo].at(e)[a] = R_NaN;
						}
						else
						{
						// Tie withdrawal contributes the opposite of tie creating
							if (pCurrentLessMissingsEtc->tieValue(e,a))
							{
								this->lstaticChangeContributions[pInfo].at(e)[a] =
									-pEffect->calculateContribution(a);
							}
							else
							{
								this->lstaticChangeContributions[pInfo].at(e)[a] =
									pEffect->calculateContribution(a);
							}
						}
					}
			}
		}
		delete pEffect;
	}
	// Restore the predictor network
	this->lpPredictorState->pNetwork(name, pPredictorNetwork);
}


/**
 * Calculates the statistics for the endowment effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkEndowmentStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	// To save a lot of unnecessary effort, check first we have some
	// endowment effects
	const vector<EffectInfo *> & rEffects =
		this->lpModel->rEndowmentEffects(pNetworkData->name());

	if (rEffects.size() > 0)
	{
		// In order to calculate the statistics of network endowment effects,
		// we need the initial network of the given period, and the network
		// of lost ties, namely, the ties that are present in the initial
		// network, not missing at the start of the period, and absent in the
		// current network.

		const Network * pInitialNetwork = pNetworkData->pNetwork(this->lperiod);
		Network * pLostTieNetwork = pInitialNetwork->clone();

		// Duplicate the current network so can overwrite the structurals
		Network * pCurrentNetwork =
			this->lpState->pNetwork(pNetworkData->name())->clone();

		// Replace values for structurally determined values from previous
		// or current period
		replaceNetwork(pCurrentNetwork,
			pNetworkData->pNetwork(this->lperiod + 1),
			pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		replaceNetwork(pCurrentNetwork,
			pNetworkData->pNetwork(this->lperiod),
			pNetworkData->pStructuralTieNetwork(this->lperiod));

		// remove missings and current

		subtractNetwork(pLostTieNetwork, pCurrentNetwork);
		subtractNetwork(pLostTieNetwork, pNetworkData->pMissingTieNetwork(this->lperiod + 1));

		// overwrite the predictor network with only start missings removed
		const Network * pPredictor =
			pNetworkData->pNetworkLessMissingStart(this->lperiod);

		// We want to pass all networks to the effects in a single state,
		// hence we overwrite the network in the predictor state.

		string name = pNetworkData->name();
		const Network * pPredictorNetwork =
			this->lpPredictorState->pNetwork(name);
		this->lpPredictorState->pNetwork(name, pPredictor);

		// Loop through the endowment effects, calculate the statistics,
		// and store them.

		EffectFactory factory(this->lpData);
		Cache cache;

		for (unsigned i = 0; i < rEffects.size(); i++)
		{
			EffectInfo * pInfo = rEffects[i];
			NetworkEffect * pEffect =
				(NetworkEffect *) factory.createEffect(pInfo);

			// Initialize the effect to work with our data and state of variables.

			pEffect->initialize(this->lpData,
				this->lpPredictorState,
				this->lperiod,
				&cache);

			if(this->lneedActorStatistics)
			{
				pair<double, double * > p = pEffect->endowmentStatistic(pLostTieNetwork, this->lneedActorStatistics);
				this->lstatistics[pInfo] = p.first;
				this->lactorStatistics[pInfo] = p.second;
			}
			else
			{
				this->lstatistics[pInfo] =
					pEffect->endowmentStatistic(pLostTieNetwork);
			}
			delete pEffect;
		}

		// Restore the predictor network
		this->lpPredictorState->pNetwork(name, pPredictorNetwork);

		delete pCurrentNetwork;
		delete pLostTieNetwork;
	}
}


/**
 * Calculates the statistics for the creation effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkCreationStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	// To save a lot of unnecessary effort, check first we have some
	// endowment effects
	const vector<EffectInfo *> & rEffects =
		this->lpModel->rCreationEffects(pNetworkData->name());

	if (rEffects.size() > 0)
	{
		// We want to pass all networks to the effects in a single state,
		// hence we overwrite the network in the predictor state.
		// We do not use the predictor network with effects this network owns.

		string name = pNetworkData->name();
		const Network * pPredictorNetwork =
			this->lpPredictorState->pNetwork(name);
		const Network * pCurrentLessMissingsEtc =
			this->lpStateLessMissingsEtc->pNetwork(name);
		this->lpPredictorState->pNetwork(name, pCurrentLessMissingsEtc);

		// In order to calculate the statistics of tie creation effects, we
		// need to know the ties that are gained, namely, those that are
		// known to be absent in the beginning of the period.

		Network * pGainedTieNetwork = pCurrentLessMissingsEtc->clone();

		subtractNetwork(pGainedTieNetwork,
			pNetworkData->pNetwork(this->lperiod));
		// not sure we need this one!
		subtractNetwork(pGainedTieNetwork,
			pNetworkData->pMissingTieNetwork(this->lperiod));

		// Loop through the creation effects, calculate the statistics,
		// and store them.

		EffectFactory factory(this->lpData);
		Cache cache;

		for (unsigned i = 0; i < rEffects.size(); i++)
		{
			EffectInfo * pInfo = rEffects[i];
			NetworkEffect * pEffect =
				(NetworkEffect *) factory.createEffect(pInfo);

			// Initialize the effect to work with our data and state of variables.

			pEffect->initialize(this->lpData,
				this->lpPredictorState,
				this->lperiod,
				&cache);

			if(this->lneedActorStatistics)
			{
				pair<double, double * > p = pEffect->creationStatistic(pGainedTieNetwork, this->lneedActorStatistics);
				this->lstatistics[pInfo] = p.first;
				this->lactorStatistics[pInfo] = p.second;
			}
			else
			{
				this->lstatistics[pInfo] =
					pEffect->creationStatistic(pGainedTieNetwork);
			}
			delete pEffect;
		}

		// Restore the predictor network
		this->lpPredictorState->pNetwork(name, pPredictorNetwork);

		delete pGainedTieNetwork;
	}
}

/**
 * Calculates the statistics for effects of the given behavior variable.
 */
void StatisticCalculator::calculateBehaviorStatistics(
	BehaviorLongitudinalData * pBehaviorData)
{
	// create a copy of the current state and zero any values missing
	// at either end of period

	const int * currentState =
		this->lpState->behaviorValues(pBehaviorData->name());

	double * currentValues  = new double[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		currentValues[i] = currentState[i] - pBehaviorData->overallMean();

		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			currentValues[i] = 0;
		}
	}

	// Construct a vector of difference values relative to previous period.
    // Values for missing values are set to 0.

	int * difference = new int[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		difference[i] =
			pBehaviorData->value(this->lperiod, i) - currentState[i];
// is this correct? centering?
		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			difference[i] = 0;
		}
	}

	// Loop through the evaluation effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEvaluationEffects =
		this->lpModel->rEvaluationEffects(pBehaviorData->name());

 	EffectFactory factory(this->lpData);
 	Cache cache;

	for (unsigned i = 0; i < rEvaluationEffects.size(); i++)
	{
		EffectInfo * pInfo = rEvaluationEffects[i];
		BehaviorEffect * pEffect =
			(BehaviorEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData,
			this->lpPredictorState,
			this->lperiod,
			&cache);

		if(this->lneedActorStatistics)
		{
			pair<double, double * > p =
				pEffect->evaluationStatistic(currentValues,this->lneedActorStatistics);
			this->lstatistics[pInfo] = p.first;
			this->lactorStatistics[pInfo] = p.second;
		}
		else
		{
			this->lstatistics[pInfo] =
				pEffect->evaluationStatistic(currentValues);
		}
		if(this->lcountStaticChangeContributions)
		{
			int  choices = 3;
			vector<double *> egosMap(pBehaviorData->n());
			this->lstaticChangeContributions.insert(make_pair(pInfo,egosMap));
			for (int e = 0; e < pBehaviorData->n(); e++)
			{
				cache.initialize(e);
				pEffect->initialize(this->lpData, this->lpPredictorState, this->lperiod, &cache);
				double * contributions = new double[choices];
				this->lstaticChangeContributions[pInfo].at(e) = contributions;
				pEffect->preprocessEgo(e);
				// no change gives no contribution
				this->lstaticChangeContributions[pInfo].at(e)[1] = 0;
				// calculate the contribution of downward change
				if ((currentState[e] > pBehaviorData->min())  // was currentValues
						&& (!pBehaviorData->upOnly(this->lperiod)))
				{
					this->lstaticChangeContributions[pInfo].at(e)[0] =
						pEffect->calculateChangeContribution(e,-1);
				}
				else
				{
					this->lstaticChangeContributions[pInfo].at(e)[0] = R_NaN;
				}
				// calculate the contribution of upward change
				if ((currentState[e] < pBehaviorData->max())  // was currentValues
						&& (!pBehaviorData->downOnly(this->lperiod)))
				{
					this->lstaticChangeContributions[pInfo].at(e)[2] =
						pEffect->calculateChangeContribution(e,1);
				}
				else
				{
					this->lstaticChangeContributions[pInfo].at(e)[2] = R_NaN;
				}
			}
		}

		delete pEffect;
	}

	// Loop through the endowment effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEndowmentEffects =
		this->lpModel->rEndowmentEffects(pBehaviorData->name());

	for (unsigned i = 0; i < rEndowmentEffects.size(); i++)
	{
		EffectInfo * pInfo = rEndowmentEffects[i];
		BehaviorEffect * pEffect =
			(BehaviorEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData,
			this->lpPredictorState,
			this->lperiod,
			&cache);

		if(this->lneedActorStatistics)
		{
			pair<double, double * > p =
				pEffect->endowmentStatistic(difference, currentValues,this->lneedActorStatistics);
			this->lstatistics[pInfo] = p.first;
			this->lactorStatistics[pInfo] = p.second;
		}
		else
		{
			this->lstatistics[pInfo] =
				pEffect->endowmentStatistic(difference, currentValues);
		}
		delete pEffect;
	}

	// Loop through the creation effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rCreationEffects =
		this->lpModel->rCreationEffects(pBehaviorData->name());

	for (unsigned i = 0; i < rCreationEffects.size(); i++)
	{
		EffectInfo * pInfo = rCreationEffects[i];
		BehaviorEffect * pEffect =
			(BehaviorEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData,
			this->lpPredictorState,
			this->lperiod,
			&cache);
		if(this->lneedActorStatistics)
		{
			pair<double, double * > p =
					pEffect->creationStatistic(difference, currentValues, this->lneedActorStatistics);
			this->lstatistics[pInfo] = p.first;
			this->lactorStatistics[pInfo] = p.second;
		}
		else
		{
			this->lstatistics[pInfo] =
				pEffect->creationStatistic(difference, currentValues);
		}
		delete pEffect;
	}

	delete[] currentValues;
	delete[] difference;
}


/**
 * Calculates the statistics for effects of the given continuous behavior
 * variable.
 */
void StatisticCalculator::calculateContinuousStatistics(
	ContinuousLongitudinalData * pContinuousData)
{
	// create a copy of the current state and zero any values missing
	// at either end of period

	const double * currentState =
		this->lpState->continuousValues(pContinuousData->name());

	double * currentValues  = new double[pContinuousData->n()];

	for (int i = 0; i < pContinuousData->n(); i++)
	{
		currentValues[i] = currentState[i]; // - pContinuousData->overallMean();

		if (pContinuousData->missing(this->lperiod, i) ||
			pContinuousData->missing(this->lperiod + 1, i))
		{
			currentValues[i] = 0;
		}
	}

	// Loop through the effects, calculate the statistics, and store them.
	const vector<EffectInfo *> & rEffects =
		this->lpModel->rEvaluationEffects(pContinuousData->name());

 	EffectFactory factory(this->lpData);
 	Cache cache;

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		ContinuousEffect * pEffect =
			(ContinuousEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData,
			this->lpPredictorState,
			this->lperiod,
			&cache);

		this->lstatistics[pInfo] =
			pEffect->evaluationStatistic(currentValues);

		delete pEffect;
	}

	delete[] currentValues;
}


/**
 * Calculates the statistics for the rate effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkRateStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	Network * pNetwork = 0;
	const Network * pConstNetwork =	this->lpStateLessMissingsEtc->
			pNetwork(pNetworkData->name());
	Network * pDifference = 0;

	// if parallel running, do processing from scratch
	if (this->lpModel->parallelRun())
	{
		// Duplicate the current network and remove those ties that are
		// missing at either end of the period. TODO set leavers back.
		// (Is the TODO not done for the current network?)

		pNetwork = this->lpState->pNetwork(pNetworkData->name())->clone();
		subtractNetwork(pNetwork,
			pNetworkData->pMissingTieNetwork(this->lperiod));
		subtractNetwork(pNetwork,
			pNetworkData->pMissingTieNetwork(this->lperiod + 1));

		// Replace values for structurally determined values from previous
		// or current period. Siena 3 does not do this ... so we won't yet

		// 	this->replaceNetwork(pNetwork,
		// 		pNetworkData->pNetwork(this->lperiod + 1),
		// 		pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		// 	this->replaceNetwork(pNetwork,
		// 		pNetworkData->pNetwork(this->lperiod),
		// 		pNetworkData->pStructuralTieNetwork(this->lperiod));

		// instead, remove all structurally determined ties at either end
		subtractNetwork(pNetwork,
			pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		subtractNetwork(pNetwork,
			pNetworkData->pStructuralTieNetwork(this->lperiod));
	}

	// construct a network of differences between current and start
	// of period.

	Network * pStart = pNetworkData->
		pNetworkLessMissing(this->lperiod)->clone();

	// must match with above code: the net result is probably equivalent
	if (this->lpModel->parallelRun())
	{
		// remove all structurally determined ties at either end
		subtractNetwork(pStart,
			pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		subtractNetwork(pStart,
			pNetworkData->pStructuralTieNetwork(this->lperiod));

		pDifference = symmetricDifference(pStart, pNetwork);
	}
	else
	{
		replaceNetwork(pStart,
			pNetworkData->pNetwork(this->lperiod + 1),
			pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		replaceNetwork(pStart,
			pNetworkData->pNetwork(this->lperiod),
			pNetworkData->pStructuralTieNetwork(this->lperiod));

		pDifference = symmetricDifference(pStart, pConstNetwork);
	}

	// basic rate distance

	if (!this->ldistances[pNetworkData])
	{
		int * array =
			new int[pNetworkData->observationCount() - 1];

		this->ldistances[pNetworkData] = array;
	}

	this->ldistances[pNetworkData][this->lperiod] = pDifference->tieCount();

	calcDifferences(pNetworkData, pDifference);

	// Loop through the rate effects, calculate the statistics,
	// and store them.
	const vector<EffectInfo *> & rEffects =
		this->lpModel->rRateEffects(pNetworkData->name());

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		//	double parameter = pInfo->parameter();
		string effectName = pInfo->effectName();
		string interactionName = pInfo->interactionName1();
		string rateType = pInfo->rateType();

		if (rateType == "covariate")
		{
			// Covariate-dependent rate effect

			//	if (parameter != 0)
			//	{
			ConstantCovariate * pConstantCovariate =
				this->lpData->pConstantCovariate(interactionName);
			ChangingCovariate * pChangingCovariate =
				this->lpData->pChangingCovariate(interactionName);
			BehaviorLongitudinalData * pBehavior =
				this->lpData->pBehaviorData(interactionName);

			if (pConstantCovariate)
			{
				double statistic = 0;

				for (TieIterator iter = pDifference->ties();
					 iter.valid();
					 iter.next())
				{
					statistic += pConstantCovariate->value(iter.ego()) *
						iter.value();
				}

				this->lstatistics[pInfo] = statistic;
			}
			else if (pChangingCovariate)
			{
				double statistic = 0;

				for (TieIterator iter = pDifference->ties();
					 iter.valid();
					 iter.next())
				{
					statistic +=
						pChangingCovariate->value(iter.ego(),
							this->lperiod) *
						iter.value();
				}

				this->lstatistics[pInfo] = statistic;
			}
			else if (pBehavior)
			{
				double statistic = 0;

				for (TieIterator iter = pDifference->ties();
					 iter.valid();
					 iter.next())
				{
					statistic +=
						pBehavior->value(this->lperiod, iter.ego()) *
						iter.value();
				}

				this->lstatistics[pInfo] = statistic;
			}
			else
			{
				throw logic_error(
					"(4) No individual covariate named '" +
					interactionName +
					"'.");
			}
			//}
		}
		else
		{
			// We expect a structural (network-dependent) rate effect here.
			NetworkLongitudinalData * pExplanatoryNetwork;
			if (interactionName == "")
			{
				pExplanatoryNetwork = pNetworkData;
			}
			else
			{
				pExplanatoryNetwork =
					this->lpData->pNetworkData(interactionName);
			}

			const Network * pStructural =
				pExplanatoryNetwork->pNetworkLessMissing(this->lperiod);

			double statistic = 0;

			for (TieIterator iter = pDifference->ties();
				 iter.valid();
				 iter.next())
			{
				if (effectName == "outRate")
				{
					statistic += pStructural->outDegree(iter.ego()) *
						iter.value();
				}
				else if (effectName == "inRate")
				{
					statistic += pStructural->inDegree(iter.ego()) *
						iter.value();
				}
				else if (effectName == "recipRate")
				{
					OneModeNetwork * pOneModeNetwork =
						(OneModeNetwork *) pStructural;
					statistic +=
						pOneModeNetwork->reciprocalDegree(iter.ego()) *
						iter.value();
				}
				else if (effectName == "outRateInv")
				{
					statistic +=
						1.0 / (pStructural->outDegree(iter.ego()) + 1) *
						iter.value();
				}
				else if (effectName == "outRateLog")
				{
					statistic +=
						log(pStructural->outDegree(iter.ego()) + 1) *
						iter.value();
				}
				else if (effectName == "inRateInv")
				{
					statistic +=
						1.0 / (pStructural->inDegree(iter.ego()) + 1) *
						iter.value();
				}
				else if (effectName == "inRateLog")
				{
					statistic +=
						log(pStructural->inDegree(iter.ego()) + 1) *
						iter.value();
				}
				else if (effectName == "recipRateInv")
				{
					OneModeNetwork * pOneModeNetwork =
						(OneModeNetwork *) pStructural;
					statistic +=
						1.0 / (pOneModeNetwork->reciprocalDegree(iter.ego()) + 1) *
						iter.value();
				}
				else if (effectName == "recipRateLog")
				{
					OneModeNetwork * pOneModeNetwork =
						(OneModeNetwork *) pStructural;
					statistic +=
						log(pOneModeNetwork->reciprocalDegree(iter.ego()) + 1) *
						iter.value();
				}
				else
				{
					throw domain_error("Unexpected rate effect " + effectName);
				}
			}

			this->lstatistics[pInfo] = statistic;
		}
	}

	delete pStart;
	delete pDifference;
	if (this->lpModel->parallelRun())
	{
		delete pNetwork;
	}

}
/**
 * Calculates the statistics for the rate effects of the given
 * behavior variable.
 */
void StatisticCalculator::calculateBehaviorRateStatistics(
	BehaviorLongitudinalData * pBehaviorData)
{
	// create a copy of the current state and zero any values missing
	//at either end of period
	const int * currentState = this->lpState->
		behaviorValues(pBehaviorData->name());

	int * currentValues  = new int[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		currentValues[i] = currentState[i];

		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			currentValues[i] = 0;
		}
	}
	// Construct a vector of absolute differences between current and start
    // of period. Differences for missing values are set to 0.

	const int * Start = pBehaviorData->values(this->lperiod);

	int * difference  = new int[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		difference[i] = abs(currentState[i] - Start[i]);
		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			difference[i] = 0;
		}
	}

	// basic rate distance

	if (!this->ldistances[pBehaviorData])
	{
		int * array =
			new int[pBehaviorData->observationCount() - 1];

		this->ldistances[pBehaviorData] = array;
	}

	int distance = 0;
	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		distance += difference[i];
	}
	this->ldistances[pBehaviorData][this->lperiod] = distance;

	// Loop through the rate effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEffects =
		this->lpModel->rRateEffects(pBehaviorData->name());

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		//	double parameter = pInfo->parameter();
		string effectName = pInfo->effectName();
		string interactionName = pInfo->interactionName1();
		string interactionName2 = pInfo->interactionName2();
		string rateType = pInfo->rateType();
		int internalEffectParameter = pInfo->internalEffectParameter();

		if (rateType == "covariate")
		{
			// Covariate-dependent rate effect

			//	if (parameter != 0)
			//	{
				ConstantCovariate * pConstantCovariate =
					this->lpData->pConstantCovariate(interactionName);
				ChangingCovariate * pChangingCovariate =
					this->lpData->pChangingCovariate(interactionName);
				BehaviorLongitudinalData * pBehavior =
					this->lpData->pBehaviorData(interactionName);
				if (pConstantCovariate)
				{
					double statistic = 0;

					for (int i = 0; i < pBehaviorData->n(); i++)
					{
						statistic += pConstantCovariate->value(i) *
							difference[i];
					}

					this->lstatistics[pInfo] = statistic;
				}
				else if (pChangingCovariate)
				{
					double statistic = 0;

					for (int i = 0; i < pBehaviorData->n(); i++)
					{
						statistic +=
							pChangingCovariate->value(i,
								this->lperiod) *
							difference[i];
					}

					this->lstatistics[pInfo] = statistic;
				}
				else if (pBehavior)
				{
					double statistic = 0;

					for (int i = 0; i < pBehaviorData->n(); i++)
					{
						statistic +=
							pBehavior->values(this->lperiod)[i] *
							difference[i];
					}

					this->lstatistics[pInfo] = statistic;
				}
				else
				{
					throw logic_error(
						"(5) No individual covariate named '" +
						interactionName +
						"'.");
				}
				//}
		}
		else if (rateType == "structural")
		{
			// We expect a structural (network-dependent) rate effect here.

			NetworkLongitudinalData *pNetworkData = this->lpData->
				pNetworkData(interactionName);
			const Network * pStructural =
				pNetworkData->pNetworkLessMissingStart(this->lperiod);

			double statistic = 0;
			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				if (effectName == "outRate")
				{
					statistic += pStructural->outDegree(i) *
						difference[i];
				}
				else if (effectName == "inRate")
				{
					statistic += pStructural->inDegree(i) *
						difference[i];
				}
				else if (effectName == "recipRate")
				{
					OneModeNetwork * pOneModeNetwork =
						(OneModeNetwork *) pStructural;
					statistic +=
						pOneModeNetwork->reciprocalDegree(i) *
						difference[i];
				}
				else if (effectName == "outRateInv")
				{
					statistic +=
						1.0 / (pStructural->outDegree(i) + 1) *
						difference[i];
				}
				else if (effectName == "outRateLog")
				{
					statistic +=
						log(pStructural->outDegree(i) + 1) *
						difference[i];
				}
				else if (effectName == "inRateInv")
				{
					statistic +=
						1.0 / (pStructural->inDegree(i) + 1) *
						difference[i];
				}
				else if (effectName == "inRateLog")
				{
					statistic +=
						log(pStructural->inDegree(i) + 1) *
						difference[i];
				}
				else if (effectName == "recipRateInv")
				{
					OneModeNetwork * pOneModeNetwork =
						(OneModeNetwork *) pStructural;
					statistic +=
						1.0 / (pOneModeNetwork->reciprocalDegree(i) + 1) *
						difference[i];
				}
				else if (effectName == "recipRateLog")
				{
					OneModeNetwork * pOneModeNetwork =
						(OneModeNetwork *) pStructural;
					statistic +=
						log(pOneModeNetwork->reciprocalDegree(i) + 1) *
						difference[i];
				}
				else
				{
					throw domain_error("Unexpected rate effect " + effectName);
				}
			}

			this->lstatistics[pInfo] = statistic;
		}
		else if (rateType == "diffusion")
		{
		    NetworkLongitudinalData *pNetworkData = this->lpData->
		        pNetworkData(interactionName);
		    const Network * pStructural =
		        pNetworkData->pNetworkLessMissingStart(this->lperiod);
			double statistic = 0;
			if (interactionName2 == "")
			{
				for (int i = 0; i < pBehaviorData->n(); i++)
				{
					if (effectName == "avExposure" ||
						effectName == "susceptAvIn" ||
						effectName == "totExposure" ||
						effectName == "infectDeg" ||
						effectName == "infectIn" ||
						effectName == "infectOut")
					{
						statistic +=
							this->calculateDiffusionRateEffect(pBehaviorData,
								pStructural, i, effectName,
								internalEffectParameter) *
													difference[i];
					}
					else
					{
						throw domain_error("Unexpected rate effect " +
							effectName);
					}
				}
			}
			else
			{
				ConstantCovariate * pConstantCovariate =
					this->lpData->pConstantCovariate(interactionName2);
				ChangingCovariate * pChangingCovariate =
					this->lpData->pChangingCovariate(interactionName2);

				for (int i = 0; i < pBehaviorData->n(); i++)
				{
					if (effectName == "susceptAvCovar" ||
						effectName == "infectCovar")
					{
						statistic +=
							this->calculateDiffusionRateEffect(pBehaviorData,
								pStructural,
								pConstantCovariate,
								pChangingCovariate, i, effectName,
								internalEffectParameter) *
													difference[i];
					}
					else
					{
						throw domain_error("Unexpected rate effect " +
							effectName);
					}
				}
			}
			this->lstatistics[pInfo] = statistic;
		}
	}

	delete[] difference;
	delete[] currentValues;
}

/**
 * Calculates the statistics for the scale effects of the given
 * continuous behavior variable.
 */
void StatisticCalculator::calculateContinuousRateStatistics(
	ContinuousLongitudinalData * pContinuousData)
{
	// create a copy of the current state and zero any values missing
	// at either end of period
	const double * currentState = this->lpState->
		continuousValues(pContinuousData->name());

	double * currentValues  = new double[pContinuousData->n()];

	for (int i = 0; i < pContinuousData->n(); i++)
	{
		currentValues[i] = currentState[i];

		if (pContinuousData->missing(this->lperiod, i) ||
			pContinuousData->missing(this->lperiod + 1, i))
		{
			currentValues[i] = 0;
		}
	}
	// Construct a vector of squared differences between current and
	// start of period. Differences for missing values are set to 0.
	const double * start = pContinuousData->values(this->lperiod);

	double * difference  = new double[pContinuousData->n()];

	for (int i = 0; i < pContinuousData->n(); i++)
	{
		difference[i] = currentState[i] - start[i];
		difference[i] *= difference[i];
		if (pContinuousData->missing(this->lperiod, i) ||
			pContinuousData->missing(this->lperiod + 1, i))
		{
			difference[i] = 0;
		}
	}

	// basic rate distance (used for estimating the SDE scale parameters)
	if (!this->lcontinuousDistances[pContinuousData])
	{
		double * array =
			new double[pContinuousData->observationCount() - 1];

		this->lcontinuousDistances[pContinuousData] = array;
	}

	double distance = 0;
	for (int i = 0; i < pContinuousData->n(); i++)
	{
		distance += difference[i];
	}
	this->lcontinuousDistances[pContinuousData][this->lperiod] = distance;

	delete[] difference;
	delete[] currentValues;
}


/**
 * Calculates the value of the diffusion rate effect for the given actor.
 */
 // note: calculateDiffusionRateEffect is also a function in DependentVariable
 // (almost the same...)
double StatisticCalculator::calculateDiffusionRateEffect(
	BehaviorLongitudinalData * pBehaviorData, const Network * pStructural,
	int i, string effectName, int internalEffectParameter)
{
	double totalAlterValue = 0;
	int numInfectedAlter = 0;
	double response = 1;
	if (pStructural->outDegree(i) > 0)
	{
		if (effectName == "avExposure")
		{
			response /= double(pStructural->outDegree(i));
		}
		else if (effectName == "susceptAvIn")
		{
			response = double(pStructural->inDegree(i)) /
				double(pStructural->outDegree(i));
		}

		for (IncidentTieIterator iter = pStructural->outTies(i);
			 iter.valid();
			 iter.next())
		{
			double alterValue = pBehaviorData->
				value(this->lperiod,iter.actor());  // this is the value at the start of the period

			if (alterValue >= 0.5)
			{
				numInfectedAlter++;
			}

			if (effectName == "infectIn")
			{
				alterValue *= pStructural->inDegree(iter.actor());
			}
			else if ((effectName == "infectOut") || (effectName == "infectDeg"))
			{
				alterValue *= pStructural->outDegree(iter.actor());
			}

			totalAlterValue += alterValue;
		}

		if (internalEffectParameter != 0)
		{
			if (numInfectedAlter < std::abs(internalEffectParameter))
			{
				totalAlterValue = 0;
			}
			else if (internalEffectParameter < 0)
			{
				if (totalAlterValue + internalEffectParameter > 0)
				{
					totalAlterValue = - internalEffectParameter;
				}
			}
		}
		totalAlterValue *= response;
	}
	return totalAlterValue;
}
/**
 * Calculates the value of the covariate dependent diffusion rate effect for
 * the given actor.
 */
double StatisticCalculator::calculateDiffusionRateEffect(
	BehaviorLongitudinalData * pBehaviorData, const Network * pStructural,
	const ConstantCovariate * pConstantCovariate,
	const ChangingCovariate * pChangingCovariate,
	int i, string effectName, int internalEffectParameter)
{
	double totalAlterValue = 0;
	double response = 1;
	int numInfectedAlter = 0;

	if (pStructural->outDegree(i) > 0)
	{
		if (effectName == "susceptAvCovar")
		{
			if (pConstantCovariate)
			{
				response = pConstantCovariate->value(i);
			}
			else if (pChangingCovariate)
			{
				response = pChangingCovariate->value(i, this->lperiod);
			}
			else
			{
				throw logic_error(
					"No individual covariate.");
			}
			response /= double(pStructural->outDegree(i));
		}

		for (IncidentTieIterator iter = pStructural->outTies(i);
			 iter.valid();
			 iter.next())
		{
			double alterValue = pBehaviorData->
				value(this->lperiod,iter.actor());

			if (alterValue >= 0.5)
			{
				numInfectedAlter++;
			}

			if (effectName == "infectCovar")
			{
				if (pConstantCovariate)
				{
					alterValue *= pConstantCovariate->value(iter.actor());
				}
				else if (pChangingCovariate)
				{
					alterValue *= pChangingCovariate->value(iter.actor(),
						this->lperiod);
				}
				else
				{
					throw logic_error("No individual covariate.");
				}
			}
			totalAlterValue += alterValue;
		}

		if (internalEffectParameter != 0)
		{
			if (numInfectedAlter < std::abs(internalEffectParameter))
			{
				totalAlterValue = 0;
			}
			else if (internalEffectParameter < 0)
			{
				if (totalAlterValue + internalEffectParameter > 0)
				{
					totalAlterValue = - internalEffectParameter;
				}
			}
		}
		totalAlterValue *= response;
	}
	return totalAlterValue;
}

void StatisticCalculator::calcDifferences(
		NetworkLongitudinalData * const pNetworkData,
		const Network* const pDifference)
{
	// names of the settings
	const std::vector<SettingInfo> & rSettingNames =
			pNetworkData->rSettingNames();
	// if there is no setting we have already stored the difference
	if (rSettingNames.empty())
	{
		return;
	}

	// create a map that stores for each edge in the difference network
	// the settings index that can explain this edge (explain means
	// the settings has the possibility to change this edge)
	std::map<pair<int, int>, std::vector<int> >* map = new std::map<
			pair<int, int>, std::vector<int> >();
	// since setting i=0 has to be the universal setting and this setting
	// can explain every edge we just store it in the map
	// TODO: This is a dependency we want to get rid off (currently user has
	// to ensure that the first setting is the universal setting ... otherwise
	// this will cause missbehaviour)
	for (TieIterator iter = pDifference->ties(); iter.valid(); iter.next())
	{
		// vector<int>(1) creates an entry with 0
		map->insert(
				make_pair(make_pair(iter.ego(), iter.alter()),
						std::vector<int>(1)));
	}
	// for each other setting (primary and covariate settings) check
	// which edges can be explained by the setting and add it to the
	// proper entry in the map
	for (unsigned i = 1; i < rSettingNames.size(); i++)
	{
		// primary
		if (i == 1)
		{
			Network * settingNetwork = pNetworkData->pNetworkLessMissing(
					this->lperiod)->clone();
			ConstantDyadicCovariate * pConstantDyadicCovariate = 0;
			ChangingDyadicCovariate * pChangingDyadicCovariate = 0;
			if (!rSettingNames[i].getCovarName().empty())
			{
				pConstantDyadicCovariate =
						this->lpData->pConstantDyadicCovariate(
								rSettingNames[i].getCovarName());
				pChangingDyadicCovariate =
						this->lpData->pChangingDyadicCovariate(
								rSettingNames[i].getCovarName());
			}

			// create a network representing primary settings
			DistanceTwoLayer primSetting;
			primSetting.onInitializationEvent(*settingNetwork);

			for (int ego = 0; ego < settingNetwork->n(); ego++)
			{
				for (IncidentTieIterator iter = settingNetwork->outTies(ego);
						iter.valid(); iter.next())
				{
					settingNetwork->setTieValue(iter.actor(), ego, 1);
				}
				for (IncidentTieIterator iter =
						primSetting.getDistanceTwoNeighbors(ego); iter.valid();
						iter.next())
				{
					settingNetwork->setTieValue(ego, iter.actor(), 1);
				}
			}
			primSetting.onNetworkDisposeEvent(*settingNetwork);
			// iterate of ties in the difference network
			for (TieIterator iter = pDifference->ties(); iter.valid();
					iter.next())
			{
				{
					// if the primary setting has this tie we know that it is
					// not existent in the simulated network and therefore we
					// can explain this edge via the primary setting
					if (settingNetwork->tieValue(iter.ego(), iter.alter())
							|| (pConstantDyadicCovariate != 0
									&& pConstantDyadicCovariate->value(
											iter.ego(), iter.alter()))
							|| (pChangingDyadicCovariate != 0
									&& pChangingDyadicCovariate->value(
											iter.ego(), iter.alter(),
											this->lperiod)))
					{
						map->find(make_pair(iter.ego(), iter.alter()))->second.push_back(
								i);
					}
				}
			}
			delete settingNetwork;
		}
		// for each covariate-defined setting, calculate
		// setting*difference network and sum.
		if (i > 1)
		{
			ConstantDyadicCovariate * pConstantDyadicCovariate =
					this->lpData->pConstantDyadicCovariate(
							rSettingNames[i].getCovarName());
			ChangingDyadicCovariate * pChangingDyadicCovariate =
					this->lpData->pChangingDyadicCovariate(
							rSettingNames[i].getCovarName());

			// iterate over the ties of the difference network
			for (TieIterator iter = pDifference->ties(); iter.valid();
					iter.next())
			{
				// if our dyadic covariate has an entry for this edge
				// we know it is not part of the simulated network and
				// we can explain it here
				if (pConstantDyadicCovariate)
				{
					if (pConstantDyadicCovariate->value(iter.ego(),
							iter.alter()))
					{
						map->find(make_pair(iter.ego(), iter.alter()))->second.push_back(
								i);
					}
				} else if (pChangingDyadicCovariate)
				{
					if (pChangingDyadicCovariate->value(iter.ego(),
							iter.alter(), this->lperiod))
					{
						map->find(make_pair(iter.ego(), iter.alter()))->second.push_back(
								i);
					}
				}
				else
				{
					throw logic_error(
							"No dyadic covariate named '"
									+ rSettingNames[i].getId() + "'.");
				}
			}
		}
	}
	// init distances with 0s
	double* distances = new double[rSettingNames.size()]();
	// init stores
	for (unsigned i = 0; i < rSettingNames.size(); i++)
	{
		if (!this->lsettingDistances[pNetworkData][rSettingNames[i].getId()])
		{
			int * array = new int[pNetworkData->observationCount() - 1];
			this->lsettingDistances[pNetworkData][rSettingNames[i].getId()] =
					array;
		}
	}
	// iterate over each entry in the map and add the distances up
	for (std::map<pair<int, int>, vector<int> >::const_iterator iter =
			map->begin(); iter != map->end(); ++iter)
	{
		// vector storing the settings that are able to explain the current tie flip
		const vector<int>& vec = iter->second;
		// the number of settings that can explain this tie flip
		double size = 1;
		if (this->lpModel->normalizeSettingRates())
		{
			size = static_cast<double>(vec.size());
		}
		// iterate over the settings and increase the distances
		for (std::vector<int>::const_iterator settingsIter = vec.begin();
				settingsIter != vec.end(); ++settingsIter)
		{
			// increase distance by 1 divided by the number of settings that
			// are able to explain this tie flip
			distances[*settingsIter] += 1 / size;
		}
	}
	// delete the map
	delete map;
	// store distances
	for (unsigned i = 0; i < rSettingNames.size(); i++)
	{
		this->lsettingDistances[pNetworkData][rSettingNames[i].getId()][this->lperiod] =
				(int) distances[i];
	}
	delete[] distances;
}

}
