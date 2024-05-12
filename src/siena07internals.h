/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07internals.h
 *
 * Description: This file contains prototypes for the internal routines
 * used to setup the data in C and create initial ML chains.
 *****************************************************************************/

#ifndef SIENA07INTERNALS_H_
#define SIENA07INTERNALS_H_

#include <Rinternals.h>

#include <vector>

namespace siena
{
	class Data;
	class Model;
	class StatisticCalculator;
	class EpochSimulation;
	class MLSimulation;
	class NetworkLongitudinalData;
	class OneModeNetworkLongitudinalData;
	class BehaviorLongitudinalData;
	class ContinuousLongitudinalData;
	class ConstantCovariate;
	class ChangingCovariate;
	class ConstantDyadicCovariate;
	class ChangingDyadicCovariate;
}

using namespace siena;

/**
 * Matches column names with indices in the effects object.
 */
void getColNos(SEXP Names, int * netTypeCol, int * nameCol, int * effectCol,
	int * parmCol, int * int1Col, int * int2Col, int * initValCol,
	int * typeCol, int * groupCol, int * periodCol, int * pointerCol,
	int * rateTypeCol, int * intptr1Col, int * intptr2Col, int * intptr3Col,
	int * settingCol);

/**
 *  updates the parameter values for each of the effects.
 */
void updateParameters(SEXP EFFECTSLIST, SEXP THETA, std::vector<Data *> *
	pGroupData, Model * pModel);

/**
 * Create one observation for a one mode Network: ties, missing, structural
 *
 */
void setupOneModeNetwork(SEXP ONEMODE,
	OneModeNetworkLongitudinalData * pNetworkData,
	int observation);

/**
 * Create all observations for a one mode Network
 *
 */
void setupOneModeObservations(SEXP ONEMODES,
	OneModeNetworkLongitudinalData *
	pOneModeNetworkLongitudinalData);

/**
 * Create one group of one mode Networks
 *
 */
void setupOneModeGroup(SEXP ONEMODEGROUP, Data * pData);

/**
 * Create one observation for a bipartite Network: ties, missing, structural
 *
 */
void setupBipartiteNetwork(SEXP BIPARTITE,
	NetworkLongitudinalData * pNetworkData,
	int observation);

/**
 * Create all observations for a bipartite Network
 *
 */
void setupBipartiteObservations(SEXP BIPARTITES,
	NetworkLongitudinalData *
	pNetworkLongitudinalData);

/**
 * Create one group of bipartite Networks
 *
 */
void setupBipartiteGroup(SEXP BIPARTITEGROUP, Data * pData);

/**
 * Create all observations for a behavior Network
 *
 */
void setupBehavior(SEXP BEHAVIOR, BehaviorLongitudinalData * pBehaviorData);


/**
 * Create one group of Behavior Networks
 *
 */
void setupBehaviorGroup(SEXP BEHGROUP, Data *pData);

/**
 * Create all observations for a continuous dependent variable
 *
 */
void setupContinuous(SEXP CONTINUOUS, ContinuousLongitudinalData * 
	pContinuousData);

/**
 * Create one group of Continuous dependent variables
 *
 */
void setupContinuousGroup(SEXP CONTGROUP, Data *pData);

/**
 * Create a constant covariate
 *
 */
void setupConstantCovariate(SEXP COCOVAR, ConstantCovariate *
	pConstantCovariate);

/**
 * Create one group of constant covariates
 *
 */
void setupConstantCovariateGroup(SEXP COCOVARGROUP, Data *pData);

/**
 * Create all observations for a changing covariate
 *
 */
void setupChangingCovariate(SEXP VARCOVAR,
	ChangingCovariate * pChangingCovariate);


/**
 * Create one group of changing covariates
 *
 */
void setupChangingCovariateGroup(SEXP VARCOVARGROUP, Data *pData);

/**
 * Create a constant dyadic covariate
 *
 */
void setupDyadicCovariate(SEXP DYADVAR,
	ConstantDyadicCovariate * pConstantDyadicCovariate);

/**
 * Create one group of constant dyadic covariates
 *
 */
void setupDyadicCovariateGroup(SEXP DYADVARGROUP, Data *pData);

/**
 * Unpack one set of values for a changing dyadic covariate
 *
 */
void unpackChangingDyadicPeriod(SEXP VARDYADVALS, ChangingDyadicCovariate *
	pChangingDyadicCovariate, int period);

/**
 * Create all observations for a changing dyadic covariate
 *
 */
void setupChangingDyadicObservations(SEXP VARDYAD,
	ChangingDyadicCovariate *
	pChangingDyadicCovariate);

/**
 * Create one group of changing dyadic covariates
 *
 */
void setupChangingDyadicCovariateGroup(SEXP VARDYADGROUP, Data * pData);

/**
 * Create the exogenous composition change events for one actor set within
 * one group.
 */
void setupExogenousEventSet(SEXP EXOGEVENTSET, Data *pData);

/**
 * Create one group of exogenous composition change events
 *
 */
void setupExogenousEventGroup(SEXP EXOGEVENTGROUP, Data *pData);

/**
 *  Creates all the basic effects for one network
 */
SEXP createEffects(SEXP EFFECTS, Model *pModel, std::vector<Data *> * pGroupData,
		const char *networkName, int effectCol,
		int parmCol, int int1Col, int int2Col,
		int initValCol, int typeCol, int groupCol,
		int periodCol, int rateTypeCol,
		int netTypeCol, int settingCol);

/**
 *  Creates all the interaction effects for one network
 */
SEXP createInteractionEffects(SEXP EFFECTS, Model *pModel,
		const char *networkName, int effectCol, int initValCol,
		int typeCol, int intptr1Col, int intptr2Col, int intptr3Col);

/**
 *  Retrieves the contributions to all possible tie flips or behavior changes
 *  for each of the effects, for one period. The call will relate to one group
 *  only, although all effects are the same apart from the basic rates. Not
 *  used in maximum likelihood.
 */
void getChangeContributionStatistics(SEXP EFFECTSLIST,
	const StatisticCalculator * pCalculator,
	std::vector<std::vector<double * > > *rChangeContributions);

/**
 *  Retrieves the statistics of individual actors for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getActorStatistics(SEXP EFFECTSLIST,
	const StatisticCalculator * pCalculator,
	std::vector<double *> *rActorStatistics);

/**
 *  Retrieves the values of the statistics and scores for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getStatistics(SEXP EFFECTSLIST,
	const StatisticCalculator * pCalculator,
	int period, int group, const Data *pData,
	const EpochSimulation * pEpochSimulation,
	std::vector<double> * rfra, std::vector<double> *rscore);

/**
 *  retrieves the values of the scores and derivatives for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Only used in maximum likelihood.
 */
void getScores(SEXP EFFECTSLIST, int period, int group,
	const MLSimulation * pMLSimulation,
	std::vector<double> * rderiv, std::vector<double> *rscore);

#endif /*SIENA07INTERNALS_H_*/
