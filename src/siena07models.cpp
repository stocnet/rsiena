/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07models.cpp
 *
 * Description: This module contains routines to interface with R,
 * running simulation models. Routines in this file are
 * visible from R.
 *****************************************************************************/
/**
 * @file
 * Runs simulations.
 */

#include <stdexcept>
#include <vector>
#include <cstring>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include "siena07internals.h"
#include "siena07utilities.h"
#include "data/Data.h"
#include "data/LongitudinalData.h"
#include "model/EffectInfo.h"
#include "model/Model.h"
#include "model/State.h"
#include "model/StatisticCalculator.h"
#include "utils/Random.h"
#include "model/EpochSimulation.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/MLSimulation.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/BehaviorChange.h"

using namespace std;
using namespace siena;

/**
 * Convert an inter SEXP to c int. If the R value isNull, return the default value.
 */
static int sexp_to_int(SEXP value, int def) {
	if (!isNull(value))
	{
		return asInteger(value);
	}
	return def;
}

extern "C"
{

/**
 *  Does one forward simulation for all the data by period within group
 */
SEXP forwardModel(SEXP DERIV, SEXP DATAPTR, SEXP SEEDS,
	SEXP FROMFINITEDIFF, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP RANDOMSEED2, SEXP RETURNDEPS, SEXP NEEDSEEDS,
	SEXP USESTREAMS, SEXP ADDCHAINTOSTORE, SEXP RETURNCHAINS, SEXP RETURNLOGLIK,
	SEXP RETURNACTORSTATISTICS, SEXP RETURNCHANGECONTRIBUTIONS)
{
	SEXP NEWRANDOMSEED = PROTECT(duplicate(RANDOMSEED2)); // for parallel testing only

	/* create a simulation and return the observed statistics and scores */

	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *) R_ExternalPtrAddr(DATAPTR);

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);
	//	Rprintf("%x %x\n", pGroupData, pModel);
	int nGroups = pGroupData->size();

	/* find total number of periods to process */
	int totObservations = totalPeriods(*pGroupData);
	int fromFiniteDiff = asInteger(FROMFINITEDIFF);
	int useStreams = asInteger(USESTREAMS);
	int addChainToStore = sexp_to_int(ADDCHAINTOSTORE, 0);
	int returnDependents = asInteger(RETURNDEPS);
	int returnChains = sexp_to_int(RETURNCHAINS, 0);
	int returnLoglik = sexp_to_int(RETURNLOGLIK, 0);
	int returnActorStatistics = sexp_to_int(RETURNACTORSTATISTICS, 0);
	int returnChangeContributions = sexp_to_int(RETURNCHANGECONTRIBUTIONS, 0);
	int deriv = asInteger(DERIV);
	int needSeeds = asInteger(NEEDSEEDS);

	/* set the deriv flag on the model */
	pModel->needScores(deriv);

	/* set the chain flag on the model */
	pModel->needChain(returnChains == 1 || addChainToStore == 1
			|| returnLoglik == 1 || returnChangeContributions == 1);

	/* set the change contribution flag on the model */
	pModel->needChangeContributions(returnChangeContributions);

	/* update the parameters */
	updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);

	/* ans will be the return value */
	SEXP ans = PROTECT(allocVector(VECSXP, 10));

	/* count up the total number of parameters */
	int dim = 0;
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		dim += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}

	/* get the random seed from R into memory */
	GetRNGstate();

	/* fra will contain the simulated statistics and must be initialised
	   to 0. Use rfra to reduce function evaluations. */
	SEXP fra;
	double * rfra;
	PROTECT(fra = allocMatrix(REALSXP, dim, totObservations));
	rfra = REAL(fra);
	for (int i = 0; i < length(fra); i++)
	{
		rfra[i] = 0;
	}
	/* ntim is the total time taken in each period (relevant for conditional
		 estimation) */
	SEXP ntim;
	double * rntim;
	PROTECT(ntim = allocVector(REALSXP, totObservations));
	rntim = REAL(ntim);
	for (int i = 0; i < length(ntim); i++)
		rntim[i] = 0.0;

	/* sims will be the returned simulated dependent variables */
	SEXP sims;
	PROTECT(sims = allocVector(VECSXP, nGroups));
	if (returnDependents)
	{
		int nVariables = (*pGroupData)[0]->rDependentVariableData().size();
		for (int group = 0; group < nGroups; group++)
		{
			SET_VECTOR_ELT(sims, group, allocVector(VECSXP, nVariables));
			for (int variable = 0; variable < nVariables; variable++)
			{
				SET_VECTOR_ELT(VECTOR_ELT(sims, group), variable,
					allocVector(VECSXP, (*pGroupData)[group]->
						observationCount() - 1));
			}
		}
	}
	/* chains will be the returned chains */
	SEXP chains = PROTECT(allocVector(VECSXP, nGroups));
	if (returnChains)
	{
		for (int group = 0; group < nGroups; group++)
		{
			SET_VECTOR_ELT(chains, group,
				allocVector(VECSXP, (*pGroupData)[group]->
					observationCount() - 1));
		}
	}

	/* changeContributionChain will be the returned change contributions for all chains */
	SEXP changeContributionChains = PROTECT(allocVector(VECSXP, nGroups));
	if (returnChangeContributions)
	{
		for (int group = 0; group < nGroups; group++)
		{
			SET_VECTOR_ELT(changeContributionChains, group,
					allocVector(VECSXP, (*pGroupData)[group]->observationCount() - 1));
		}
	}

	/* actorStats will be the returned statistics of individual actors*/
	SEXP actorStats = PROTECT(allocVector(VECSXP,nGroups));
	if(returnActorStatistics)
	{
		SEXP NETWORKTYPES;
		NETWORKTYPES = createRObjectAttributes(EFFECTSLIST, actorStats);
		int objEffects = length(NETWORKTYPES);
		for (int group = 0; group < nGroups; group++)
		{
			SET_VECTOR_ELT(actorStats, group, allocVector(VECSXP, (*pGroupData)[group]->observationCount()-1));
			for (int p = 0; p < (*pGroupData)[group]->observationCount()-1; p++)
			{
				SET_VECTOR_ELT(VECTOR_ELT(actorStats,group), p, allocVector(VECSXP,objEffects));
			}
		}
	}

	/* loglik will be the returned log likelihoods */
	SEXP logliks = PROTECT(allocVector(REALSXP, totObservations));
	double * Rlogliks = REAL(logliks);

	/* seed store is a list to save the random states */
	SEXP seedstore = PROTECT(allocVector(VECSXP, nGroups));
	for (int group = 0; group < nGroups; group++)
	{
		SET_VECTOR_ELT(seedstore, group,
			allocVector(VECSXP, (*pGroupData)[group]->observationCount() - 1));
	}

	/* rs will allow us to access or set the .Random.seed in R */
	SEXP rs = PROTECT(install(".Random.seed"));

	/* scores will hold the return values of the scores */
	SEXP scores;
	double *rscores;
	PROTECT(scores = allocMatrix(REALSXP, dim, totObservations));
	rscores = REAL(scores);
	for (int i = 0; i < length(scores); i++)
		rscores[i] = 0.0;

	int periodFromStart = 0;

	SEXP Cgstr = R_NilValue;
	SEXP STREAMS = R_NilValue;
	SEXP ans2, ans3, ans4, R_fcall1, R_fcall2, R_fcall3, R_fcall4;
	SEXP seedvector;

	if (useStreams)
	{
		// create an R character string
		PROTECT(Cgstr = allocVector(STRSXP,1));
		SET_STRING_ELT(Cgstr, 0, mkChar("Cg"));

		// find out which stream we are using
		PROTECT(R_fcall1 = lang1(install(".lec.GetStreams")));
		PROTECT(STREAMS = eval(R_fcall1, R_GlobalEnv));
	}

	/* group loop here */
	for (int group = 0; group < nGroups; group++)
	{
		/* random states need store (not fromFiniteDiff)
		   or restore (fromFiniteDiff) for each period
		   within each  group */
		SEXP seeds = R_NilValue;
		if (fromFiniteDiff)
		{
			seeds = VECTOR_ELT(SEEDS, group);
		}

		/* find out how many periods in this Data object */
		Data * pData = (*pGroupData)[group];
		int observations = pData->observationCount();

		/* create my epochsimulation object */
		EpochSimulation * pEpochSimulation  = new
			EpochSimulation(pData, pModel);

		for (int period = 0; period < observations - 1; period++)
		{

			periodFromStart++;

			if (!isNull(RANDOMSEED2)) /* parallel testing versus Siena3 */
			{
				// overwrite R's random number seed
				defineVar(rs, RANDOMSEED2, R_GlobalEnv);
				// get it into memory
				GetRNGstate();
				// move on one
				nextDouble();
				// write it back to R
				PutRNGstate();
			}
			else /* normal run */
			{
				if (fromFiniteDiff) /* restore state */
				{
					if (useStreams) /* using lecuyer random numbers */
					{
						// overwrite the current state in R
						PROTECT(R_fcall2 = lang4(install("[[<-"),
								install(".lec.Random.seed.table"), Cgstr,
								VECTOR_ELT(seeds, period)));
						PROTECT(ans2 = eval(R_fcall2, R_GlobalEnv));
						// get the overwritten state into C table
						PROTECT(R_fcall3 =
							lang2(install(".lec.CurrentStream"),
								STREAMS));
						PROTECT(ans3 = eval(R_fcall3, R_GlobalEnv));
						UNPROTECT(4);
					}
					else /* using normal random numbers */
					{
						// overwrite R's current state
						defineVar(rs, VECTOR_ELT(seeds, period),
							R_GlobalEnv);
						// get the value from .Random.seed into memory
						GetRNGstate();
					}
				}
				else /* save state */
				{
					if (needSeeds)
					{
						if (useStreams)
						{
							PROTECT(R_fcall2 =
								lang2(install(".lec.ResetNextSubstream"),
									STREAMS));
							PROTECT(ans2 = eval(R_fcall2, R_GlobalEnv));

							PROTECT(R_fcall3 =
								lang2(install(".lec.CurrentStream"),
									STREAMS));
							PROTECT(ans3 = eval(R_fcall3, R_GlobalEnv));
							// get the relevant current state from R
							PROTECT(R_fcall4 = lang3(install("[["),
									install(".lec.Random.seed.table"),
									Cgstr));
							ans4 = eval(R_fcall4, R_GlobalEnv);
							// value is not kept unless we duplicate it
							PROTECT(seedvector = duplicate(ans4));
							// store the Cg values
							SET_VECTOR_ELT(VECTOR_ELT(seedstore, group),
								period, seedvector);
							UNPROTECT(6);
						}
						else
						{
							PutRNGstate();
							SET_VECTOR_ELT(VECTOR_ELT(seedstore, group),
								period, findVar(rs, R_GlobalEnv));
						}
					}
				}
			}
			/* only needed for forward chains */
			if (pModel->needChain())
			{
				pEpochSimulation->clearChain();
			}
			/* run the epoch simulation for this period */
			pEpochSimulation->runEpoch(period);
			State State(pEpochSimulation);
			StatisticCalculator Calculator(pData, pModel, &State,
				period);
			vector<double> statistic(dim);
			vector<double> score(dim);
			getStatistics(EFFECTSLIST, &Calculator,
				period, group, pData, pEpochSimulation,
				&statistic, &score);  // ABC
			/* fill up matrices for  return value list */
			int iii = (periodFromStart - 1) * dim;
			for (unsigned effectNo = 0; effectNo < statistic.size();
				 effectNo++)
			{
				rfra[iii + effectNo] = statistic[effectNo];
				rscores[iii + effectNo] = score[effectNo];
			}
			if(returnActorStatistics)
			{
				StatisticCalculator Calculator(pData, pModel, &State, period, returnActorStatistics);
				vector<double *> actorStatistics;
				getActorStatistics(EFFECTSLIST, &Calculator, &actorStatistics);
				int actors = pData->rDependentVariableData()[0]->n();
				for(unsigned e = 0; e < actorStatistics.size(); e++)
				{
					SEXP actorStatsValues;
					PROTECT(actorStatsValues = allocVector(REALSXP, actors));
					double * astats = REAL(actorStatsValues);
					for(int i = 0; i < actors; i++)
					{
						astats[i] = actorStatistics.at(e)[i];
					}
					SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(actorStats,group), period), e, actorStatsValues);
					UNPROTECT(1);
				}
			}
			if (pModel->conditional())
			{
				rntim[periodFromStart - 1] = pEpochSimulation->time();
			}
			// get simulated network
			if (returnDependents)
			{
				const vector<DependentVariable *> rVariables =
					pEpochSimulation->rVariables();
				for (unsigned i = 0; i < rVariables.size(); i++)
				{
					// Rprintf("attach var %d (2)\n", i);
					SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(sims, group), i), period,
							var_to_sexp(rVariables[i]));
				}
			}
			if (returnChains)
			{

				SEXP thisChain =
					getChainList(*(pEpochSimulation->pChain()));

				SET_VECTOR_ELT(VECTOR_ELT(chains, group), period,
					thisChain);

			}
			if(returnChangeContributions)
			{
				SEXP thisChangeContributions = getChangeContributionsList(*(pEpochSimulation->pChain()), EFFECTSLIST);
				SET_VECTOR_ELT(VECTOR_ELT(changeContributionChains, group), period, thisChangeContributions);
			}
			if (returnLoglik)
			{
				pEpochSimulation->simpleRates(pModel->simpleRates());
				Rlogliks[periodFromStart - 1] =
					pEpochSimulation->calculateLikelihood();
			}
			if (addChainToStore)
			{
				pModel->chainStore(*(pEpochSimulation->pChain()),
					periodFromStart - 1);
			}
		} /* end of period */
		delete pEpochSimulation;
	} /* end of group */

	/* send the .Random.seed back to R */
	PutRNGstate();
	NEWRANDOMSEED = findVar(rs, R_GlobalEnv);

	/* set up the return object */
	if (!fromFiniteDiff)
	{
		if (needSeeds)
		{
			SET_VECTOR_ELT(ans, 2, seedstore);
		}
	}
	if (deriv)
	{
		SET_VECTOR_ELT(ans, 1, scores);
	}
	if (returnDependents)
	{
		SET_VECTOR_ELT(ans, 5, sims);/* not done in phase 2 !!test this */
	}
	SET_VECTOR_ELT(ans, 0, fra);
	SET_VECTOR_ELT(ans, 3, ntim);

	if (!isNull(RANDOMSEED2))
	{
		SET_VECTOR_ELT(ans, 4, NEWRANDOMSEED);
	}
	if (useStreams)
	{
		UNPROTECT(3);
	}
	SET_VECTOR_ELT(ans, 6, chains);
	SET_VECTOR_ELT(ans, 7, logliks);
	SET_VECTOR_ELT(ans, 8, changeContributionChains);
	SET_VECTOR_ELT(ans, 9, actorStats);
	UNPROTECT(12);
	return(ans);
}

/** Does some MH steps for a specified group and period.
 * Designed to be used for parallel processing, and currently the only
 * function available. Loop is always constructed in R. Probably would be
 * better to do it in C unless parallel processing.
 */
SEXP mlPeriod(SEXP DERIV, SEXP DATAPTR, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP GROUP, SEXP PERIOD,
	SEXP NRUNMH, SEXP ADDCHAINTOSTORE,
	SEXP RETURNDATAFRAME, SEXP RETURNDEPS, SEXP RETURNCHAINS,
	SEXP RETURNLOGLIK, SEXP ONLYLOGLIK)
{
	/* do some MH steps and return some or all of the scores and derivs
	   of the chain at the end or the calculated log likelihood*/

	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);

	int group = asInteger(GROUP) - 1;
	int period = asInteger(PERIOD) - 1;
	int groupPeriod = periodFromStart(*pGroupData, group, period);

	/* get hold of the data object */
	Data * pData = (*pGroupData)[group];

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	/* update the parameters */
	updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);

	/* create the ML simulation object */
	MLSimulation * pMLSimulation = new MLSimulation(pData, pModel);

	/* initialize some data from model */
	pMLSimulation->simpleRates(pModel->simpleRates());

	pMLSimulation->currentPermutationLength(
		pModel->currentPermutationLength(period));

	// next calls are ambiguous unless I use a const pModel
	const Model * pConstModel = pModel;

	pMLSimulation->
		missingNetworkProbability(pConstModel->
			missingNetworkProbability(groupPeriod));
	pMLSimulation->
		missingBehaviorProbability(pConstModel->
			missingBehaviorProbability(groupPeriod));

	// get chain for this period from model
	Chain * pChain = pModel->rChainStore(groupPeriod).back();

	// then copy the chain to the MLSimulation object. (deleting new one first)
	pMLSimulation->pChain(pChain->copyChain());
	//	Rprintf(" %d\n", pMLSimulation->pChain()->ministepCount());

	int addChainToStore = 0;
	if (!isNull(ADDCHAINTOSTORE))
	{
		addChainToStore = asInteger(ADDCHAINTOSTORE);
	}
	// prepare to recreate after the simulation
	if (!addChainToStore)
	{
		pModel->deleteLastChainStore(groupPeriod);
	}

	int returnDeps = 0;
	if (!isNull(RETURNDEPS))
	{
		returnDeps = asInteger(RETURNDEPS);
	}
	int returnChains = 0;
	if (!isNull(RETURNCHAINS))
	{
		returnChains = asInteger(RETURNCHAINS);
	}
	int returnDataFrame = 0;
	if (!isNull(RETURNDATAFRAME))
	{
		returnDataFrame = asInteger(RETURNDATAFRAME);
	}

	int deriv = asInteger(DERIV);
	int returnLoglik = 0;
	if (!isNull(RETURNLOGLIK))
	{
		returnLoglik = asInteger(RETURNLOGLIK);
	}
	int onlyLoglik = 0;
	if (!isNull(ONLYLOGLIK))
	{
		onlyLoglik = asInteger(ONLYLOGLIK);
	}

	// get the value from .Random.seed into memory
	GetRNGstate();

	pModel->needScores(false);
	pModel->needDerivatives(false);

	int nrunMH = asInteger(NRUNMH);
	pModel->numberMLSteps(nrunMH);

	/* run the epoch simulation for this period */
	pMLSimulation->runEpoch(period);

	/* run through current state of chain and calculate
	   scores and derivatives */
	pModel->needScores(!onlyLoglik);
	pModel->needDerivatives(deriv);

	pMLSimulation->updateProbabilities(pMLSimulation->pChain(),
		pMLSimulation->pChain()->pFirst()->pNext(),
		pMLSimulation->pChain()->pLast()->pPrevious());

	double loglik = 0;
	if (returnLoglik)
	{
		loglik =
			pMLSimulation->calculateLikelihood();
	}

 	/* store chain on Model */
	pChain = pMLSimulation->pChain();
	pChain->createInitialStateDifferences();
	pMLSimulation->createEndStateDifferences();
	pModel->chainStore(*pChain, groupPeriod);

	/* and current permutation length */
	pModel->currentPermutationLength(period,
		pMLSimulation->currentPermutationLength());

	/* get hold of the statistics for accept and reject */
	const vector < DependentVariable * > & rVariables =
		pMLSimulation->rVariables();
	int numberVariables = rVariables.size();

	SEXP accepts;
	PROTECT(accepts = allocMatrix(INTSXP, numberVariables, 9));
	SEXP rejects;
	PROTECT(rejects = allocMatrix(INTSXP, numberVariables, 9));
	SEXP aborts;
	PROTECT(aborts = allocVector(INTSXP, 9));
	int * iaccepts = INTEGER(accepts);
	int * irejects = INTEGER(rejects);
	int * iaborts = INTEGER(aborts);
	for (int i = 0; i < 9; i++)
	{
		iaborts[i] = pMLSimulation->aborted(i);
		for (int j = 0; j < numberVariables; j++)
		{
			iaccepts[i + 9 * j] = rVariables[j]->acceptances(i);
			irejects[i + 9 * j] = rVariables[j]->rejections(i);
		}
	}

	/* sims will be the returned chain */
	SEXP sims = PROTECT(allocVector(VECSXP, 1));
	/* theseDeps will be the returned dependent variables */
	SEXP theseDeps = PROTECT(allocVector(VECSXP, numberVariables));
	int nProtects = 5;

	if (returnDeps)
	{
		// get simulated last state
		pMLSimulation->gotoLastState();
		const vector<DependentVariable *> rVariables = pMLSimulation->rVariables();
		for (unsigned i = 0; i < numberVariables; i++)
		{
			SET_VECTOR_ELT(theseDeps, i, var_to_sexp(rVariables[i]));
		}
	}

	// get chain
	if (returnChains)
	{
		SEXP theseValues;
		if (returnDataFrame)
		{
			PROTECT(theseValues =
				duplicate(getChainDFPlus(*(pMLSimulation->pChain()), true)));
		}
		else
		{
			PROTECT(theseValues =
				duplicate(getChainList(*(pMLSimulation->pChain()))));
		}
		nProtects++;
		SET_VECTOR_ELT(sims, 0, theseValues);
	}

	/* ans will be the return value */
	SEXP ans;

	if (!onlyLoglik)
	{
		/* count up the total number of parameters */
		int dim = 0;
		for (int i = 0; i < length(EFFECTSLIST); i++)
		{
			dim += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
		}


		/* dff will hold the return values of the derivatives */
		SEXP dff = R_NilValue;
		double *rdff = 0;
		if (deriv)
		{
			PROTECT(dff = allocVector(REALSXP, dim * dim));
			nProtects++;
			rdff = REAL(dff);
			for (int i = 0; i < length(dff); i++)
			{
				rdff[i] = 0.0;
			}
		}
		/* collect the scores and derivatives */
		vector<double> derivs(dim * dim);
		vector<double> score(dim);

		getScores(EFFECTSLIST, 	period, group, pMLSimulation,
			&derivs, &score);

		/* fra will contain the scores and must be initialised
		   to 0. Use rfra to reduce function evaluations. */
		SEXP fra;
		double * rfra;
		PROTECT(fra = allocVector(REALSXP, dim));
		nProtects++;
		rfra = REAL(fra);
		for (int i = 0; i < length(fra); i++)
		{
			rfra[i] = 0;
		}
		/* fill up vectors for  return value list */
		for (unsigned effectNo = 0; effectNo < score.size();
			 effectNo++)
		{
			rfra[effectNo] = score[effectNo];
		}
		if (deriv)
		{
			for (unsigned ii = 0; ii < derivs.size(); ii++)
			{
				rdff[ii] = derivs[ii];
			}
		}
		PROTECT(ans = allocVector(VECSXP, 12));
		nProtects++;

		/* set up the return object */
		if (deriv)
		{
			SET_VECTOR_ELT(ans, 6, dff);
		}
		if (returnChains)
		{
			SET_VECTOR_ELT(ans, 5, sims);/* not done in phase 2 !!!!test this*/
		}
		SET_VECTOR_ELT(ans, 0, fra);
		SET_VECTOR_ELT(ans, 7, accepts);
		SET_VECTOR_ELT(ans, 8, rejects);
		SET_VECTOR_ELT(ans, 9, aborts);
		SET_VECTOR_ELT(ans, 10, ScalarReal(loglik));
		SET_VECTOR_ELT(ans, 11, theseDeps);

//	PrintValue(getChainDF(*pChain, true));
	}
	else
	{
		PROTECT(ans = allocVector(VECSXP, 5));
		nProtects++;
		SET_VECTOR_ELT(ans, 0, ScalarReal(loglik));
		SET_VECTOR_ELT(ans, 1, accepts);
		SET_VECTOR_ELT(ans, 2, rejects);
		SET_VECTOR_ELT(ans, 3, aborts);
//		SET_VECTOR_ELT(ans, 4, sims);
//		SET_VECTOR_ELT(ans, 5, theseDeps);
	}

	PutRNGstate();
	delete pMLSimulation;
	UNPROTECT(nProtects);
	return(ans);
}

/**
 * Clears the chains that have been stored on a model for a particular period
 * from start,  Leave KEEP ones for reuse
 */
SEXP clearStoredChains(SEXP MODELPTR, SEXP KEEP, SEXP GROUPPERIOD)
{
	int keep = asInteger(KEEP);
	int groupPeriod = asInteger(GROUPPERIOD) - 1;

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);
	pModel->clearChainStore(keep, groupPeriod);

	return R_NilValue;
}



}
