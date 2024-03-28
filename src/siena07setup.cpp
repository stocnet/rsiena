/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07setup.cpp
 *
 * Description: This module contains routines to interface with R,
 * setting up the Data object with data from R. Routines in this file are
 * visible from R.
 *****************************************************************************/
/**
 * @file
 * Sets up the Data object with data from R
 */


#include <vector>
#include <cstring>
#include <exception>
#include <R_ext/Random.h>
#include <R_ext/Print.h>

#include "siena07setup.h"
#include "siena07internals.h"
#include "siena07utilities.h"
#include "data/Data.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/Model.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/State.h"
#include "model/StatisticCalculator.h"
#include "model/EffectInfo.h"
#include "data/ActorSet.h"
#include "model/ml/MLSimulation.h"
#include "model/variables/DependentVariable.h"
#include <R_ext/Error.h>
#include <Rinternals.h>

using namespace std;
using namespace siena;

/**
 * Convert an inter SEXP to c int. If the R value isNull, return the default value.
 */
static int sexp_to_int(SEXP value, int def) {
	if (!Rf_isNull(value))
	{
		return Rf_asInteger(value);
	}
	return def;
}

extern "C"
{

/**
 *  Creates an array of pointers to Data objects, one for each group
 *  and returns the address of the array to R. Also creates the actor sets
 *  for each group.
 */
SEXP setupData(SEXP OBSERVATIONSLIST, SEXP ACTORSLIST)
{
	/* make error messages go back to R nicely */
	set_terminate(Rterminate);

	int nGroups = Rf_length(OBSERVATIONSLIST);

	vector<Data *> *pGroupData = new vector <Data *>;

	for (int group = 0; group < nGroups; group++)
	{
		int observations = INTEGER(VECTOR_ELT(OBSERVATIONSLIST, group))[0];

		pGroupData->push_back(new Data(observations));
		int nNodeSets = Rf_length(VECTOR_ELT(ACTORSLIST, group));

		for (int nodeSet = 0; nodeSet < nNodeSets; nodeSet++)
		{
			SEXP nsn;
			PROTECT(nsn = Rf_install("nodeSetName"));
			SEXP nodeSetName = Rf_getAttrib(VECTOR_ELT(VECTOR_ELT(ACTORSLIST,
						group),
					nodeSet), nsn);
			(*pGroupData)[group]->
				createActorSet(CHAR(STRING_ELT(nodeSetName, 0)),
					Rf_length(VECTOR_ELT(VECTOR_ELT(ACTORSLIST,
								group), nodeSet)));
			UNPROTECT(1);
		}
	}
	SEXP RpData;
	RpData = R_MakeExternalPtr((void *) pGroupData, R_NilValue,
		R_NilValue);
	return RpData;
}

/**
 *  Creates all the groups of one mode networks in the data
 *
 */
SEXP OneMode(SEXP RpData, SEXP ONEMODELIST)
{
	/* retrieve the address of our data */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);

	int nGroups = pGroupData->size();
	/* one mode networks are passed in as list of edgelists with attributes
	   giving the size of the network */
	if (nGroups != Rf_length(ONEMODELIST) )
	{
		Rf_error("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupOneModeGroup(VECTOR_ELT(ONEMODELIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}

/**
 *  Creates all the groups of bipartite networks in the data
 *
 */
SEXP Bipartite(SEXP RpData, SEXP BIPARTITELIST)
{
/* retrieve the address of our data */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* bipartite networks are passed in as list of edgelists with attributes
   giving the size of the network */
	if (nGroups != Rf_length(BIPARTITELIST) )
	{
		Rf_error("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupBipartiteGroup(VECTOR_ELT(BIPARTITELIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}

/**
 *  Creates all the groups of behavior networks in the data
 */
SEXP Behavior(SEXP RpData, SEXP BEHLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* behavior networks are passed in a list of lists of two matrices,
   one of values, one of missing values (boolean) */
	if (nGroups != Rf_length(BEHLIST) )
	{
		Rf_error("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupBehaviorGroup(VECTOR_ELT(BEHLIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}

/**
 *  Creates all the groups of continuous dependent variables in the data
 */
SEXP Continuous(SEXP RpData, SEXP CONTLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* continuous dependent variables are passed in a list of lists of two 
   matrices, one of values, one of missing values (boolean) */
	if (nGroups != Rf_length(CONTLIST) )
	{
		Rf_error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupContinuousGroup(VECTOR_ELT(CONTLIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}

/**
 *  Creates all the groups of constant covariates in the data
 */
SEXP ConstantCovariates(SEXP RpData, SEXP COCOVARLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* constant covariates are passed in as vectors with embedded missing values */
/* ignore missings for now */
	if (nGroups != Rf_length(COCOVARLIST) )
	{
		Rf_error("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupConstantCovariateGroup(VECTOR_ELT(COCOVARLIST,
				group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of changing covariates in the data
 */

SEXP ChangingCovariates(SEXP RpData, SEXP VARCOVARLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* changing covariates are passed in as matrices with embedded missing values */
/* ignore missings for now */
	if (nGroups != Rf_length(VARCOVARLIST) )
	{
		Rf_error("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupChangingCovariateGroup(VECTOR_ELT(VARCOVARLIST,
				group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of constant dyadic covariates in the data
 */

SEXP DyadicCovariates(SEXP RpData, SEXP DYADVARLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* dyadic covariates are passed in as edgelists with embedded missing values */
/* ignore missings for now */
	if (nGroups != Rf_length(DYADVARLIST) )
	{
		Rf_error("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupDyadicCovariateGroup(VECTOR_ELT(DYADVARLIST,
				group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of changing dyadic covariates in the data
 */

SEXP ChangingDyadicCovariates(SEXP RpData, SEXP VARDYADLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* dyadic covariates are passed in as lists of edgelists
   with embedded missing values */
/* ignore missings for now */
	if (nGroups != Rf_length(VARDYADLIST) )
	{
		Rf_error("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupChangingDyadicCovariateGroup(VECTOR_ELT(VARDYADLIST,
				group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the composition change events in the data
 */
SEXP ExogEvent(SEXP RpData, SEXP EXOGEVENTLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();

	/* data for each actor set in each group consists of
	   two items: first a list of events: event type, period, actor, time.
	   Secondly a matrix of booleans indicating whether active at the start
	   of the period. Final period exists in the latter but probably is not
	   necessary. */

	if (nGroups != Rf_length(EXOGEVENTLIST) )
	{
		Rf_error("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupExogenousEventGroup(VECTOR_ELT(EXOGEVENTLIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}

/**
 *  Sets the pairwise constraints for the data
 */
SEXP Constraints(SEXP RpData, SEXP FROMHIGHERLIST, SEXP TOHIGHERLIST,
	SEXP FROMDISJOINTLIST, SEXP TODISJOINTLIST,
	SEXP FROMATLEASTONELIST, SEXP TOATLEASTONELIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();

	for (int group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];

		/* higher */
		for (int i = 0; i < Rf_length(FROMHIGHERLIST); i++)
		{
			pData->
				addNetworkConstraint(CHAR(STRING_ELT(FROMHIGHERLIST, i)),
					CHAR(STRING_ELT(TOHIGHERLIST, i)), HIGHER);
		}
		/* disjoint */
		for (int i = 0; i < Rf_length(FROMDISJOINTLIST); i++)
		{
			pData->
				addNetworkConstraint(CHAR(STRING_ELT(FROMDISJOINTLIST, i)),
					CHAR(STRING_ELT(TODISJOINTLIST, i)), DISJOINT);
		}
		/* at least one */
		for (int i = 0; i < Rf_length(FROMATLEASTONELIST); i++)
		{
			pData->
				addNetworkConstraint(CHAR(STRING_ELT(FROMATLEASTONELIST,
							i)),
					CHAR(STRING_ELT(TOATLEASTONELIST, i)), AT_LEAST_ONE);
		}
	}
	return R_NilValue;
}


/**
 *  creates the requested basic effects
 */

SEXP effects(SEXP RpData, SEXP EFFECTSLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);

	Model * pModel = new Model();

	// for the SDE part of the model, in particular the scale parameters,
	// it is important that numberOfPeriods in the model is set very early 
	// on (waiting till setupModelOptions is too late); however this is 
	// still only implemented for a single group
	if(pGroupData->size() > 1)
	{
		//error("SDE not implemented yet for multiple groups (note nr. of obs!)");
		pModel->numberOfPeriods(totalPeriods(*pGroupData));
	}
	else // single group
	{
		pModel->numberOfPeriods((*(*pGroupData)[0]).observationCount() - 1);
	}
	
	/* find the number of columns of the data frame (all will be the same
	   as they are split in R just before the call) */
	//	int n = Rf_length(VECTOR_ELT(EFFECTSLIST, 0));
	// get the column names from the names attribute
	SEXP cols;
	PROTECT(cols = Rf_install("names"));
	SEXP Names = Rf_getAttrib(VECTOR_ELT(EFFECTSLIST, 0), cols);

	int netTypeCol; /* net type */
	int nameCol; /* network name */
	int effectCol;  /* short name of effect */
	int parmCol;
	int int1Col;
	int int2Col;
	int initValCol;
	int typeCol;
	int groupCol;
	int periodCol;
	int pointerCol;
	int rateTypeCol;
	int intptr1Col;
	int intptr2Col;
	int intptr3Col;
	int settingCol;

// Get the column numbers:
	getColNos(Names, &netTypeCol, &nameCol, &effectCol,
		&parmCol, &int1Col, &int2Col, &initValCol,
		&typeCol, &groupCol, &periodCol, &pointerCol,
		&rateTypeCol, &intptr1Col, &intptr2Col, &intptr3Col,
		&settingCol);

	/* create a structure for the return values */
	SEXP pointers;
	PROTECT(pointers = Rf_allocVector(VECSXP, Rf_length(EFFECTSLIST)));

	/* loop over the different dependent variables 
	 * in case there are continuous variables in the model, the 
	 * length of the effectlist is nrVar + 1 (extra var "sde") */
	for (int i = 0; i < Rf_length(EFFECTSLIST); i++)
	{
		const char * networkName =  CHAR(STRING_ELT(
				VECTOR_ELT(VECTOR_ELT(
						EFFECTSLIST, i),
					nameCol), 0));

		SEXP ptrs =
			createEffects(VECTOR_ELT(EFFECTSLIST, i), pModel, pGroupData,
				networkName, effectCol, parmCol, int1Col,
				int2Col, initValCol, typeCol, groupCol,
				periodCol, rateTypeCol, netTypeCol,
				settingCol);

		SET_VECTOR_ELT(pointers, i, ptrs);

	}
	SEXP RpModel;
	PROTECT (RpModel = Rf_allocVector(VECSXP, 1));
	SET_VECTOR_ELT(RpModel, 0, R_MakeExternalPtr((void *) pModel,
			R_NilValue,
			R_NilValue));


	/* ans will be the return value */
	SEXP ans;
	PROTECT(ans = Rf_allocVector(VECSXP, 2));

	SET_VECTOR_ELT(ans, 1, pointers);
	SET_VECTOR_ELT(ans, 0, RpModel);

	UNPROTECT(4);

	return ans;
}
/**
 *  creates the requested interaction effects
 */

SEXP interactionEffects(SEXP RpModel, SEXP EFFECTSLIST)
{
	Model * pModel = (Model *) R_ExternalPtrAddr(RpModel);

	// get the column names from the names attribute

	SEXP cols;
	PROTECT(cols = Rf_install("names"));
	SEXP Names = Rf_getAttrib(VECTOR_ELT(EFFECTSLIST, 0), cols);

	int netTypeCol; /* net type */
	int nameCol; /* network name */
	int effectCol;  /* short name of effect */
	int parmCol;
	int int1Col;
	int int2Col;
	int initValCol;
	int typeCol;
	int groupCol;
	int periodCol;
	int pointerCol;
	int rateTypeCol;
	int intptr1Col;
	int intptr2Col;
	int intptr3Col;
	int settingCol;

// Get the column numbers:
	getColNos(Names, &netTypeCol, &nameCol, &effectCol,
		&parmCol, &int1Col, &int2Col, &initValCol,
		&typeCol, &groupCol, &periodCol, &pointerCol,
		&rateTypeCol, &intptr1Col, &intptr2Col, &intptr3Col, &settingCol);

	/* create a structure for the return values */
	SEXP pointers;
	PROTECT(pointers = Rf_allocVector(VECSXP, Rf_length(EFFECTSLIST)));

	/* loop over the different dependent variables */
	for (int i = 0; i < Rf_length(EFFECTSLIST); i++)
	{
		if (Rf_length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0)) > 0)
		{
			const char * networkName =
				CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i),
							nameCol), 0));

			SEXP ptrs =
				createInteractionEffects(VECTOR_ELT(EFFECTSLIST, i),
					pModel,	networkName, effectCol, initValCol, typeCol,
					intptr1Col, intptr2Col, intptr3Col);

			SET_VECTOR_ELT(pointers, i, ptrs);
		}
		else
		{
			SET_VECTOR_ELT(pointers, i,
				R_MakeExternalPtr((void *) 0,
					R_NilValue, R_NilValue));
		}
	}
	/* ans will be the return value */
	SEXP ans;
	PROTECT(ans = Rf_allocVector(VECSXP, 1));

	SET_VECTOR_ELT(ans, 0, pointers);

	UNPROTECT(3);

	return ans;
}
/**
 *  removes the objects created for the data.
 */

SEXP deleteData(SEXP RpData)
{

	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	vector<Data *>::iterator it = pGroupData->begin();
	while(it != pGroupData->end())
	{
		delete *it;
		pGroupData->erase(it); /* not sure I need to do this */
	}
	//	Rprintf("%d delete\n", pGroupData->size());
	delete pGroupData;
	return R_NilValue;
}
/**
 *  removes the model object.
 */

SEXP deleteModel(SEXP RpModel)
{
	Model * pModel = (Model *) R_ExternalPtrAddr(RpModel);
	delete pModel;
	//	Rprintf("deleteModel\n");
	return R_NilValue;
}

/**
 *  sets up the model options of MAXDEGREE, UNIVERSALOFFSET, CONDITIONAL
 */
SEXP setupModelOptions(SEXP DATAPTR, SEXP MODELPTR, SEXP MAXDEGREE,
	SEXP UNIVERSALOFFSET,
	SEXP CONDVAR, SEXP CONDTARGETS, SEXP PROFILEDATA, SEXP PARALLELRUN,
	SEXP MODELTYPE, SEXP BEHMODELTYPE, SEXP SIMPLERATES, SEXP NORMSETRATES)
{
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);
	int nGroups = pGroupData->size();

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	if(!Rf_isNull(NORMSETRATES)){
		pModel->normalizeSettingRates(*(LOGICAL(NORMSETRATES)));
	}

	int totObservations = totalPeriods(*pGroupData);

	pModel->numberOfPeriods(totObservations);

	if (!Rf_isNull(CONDVAR))
	{
		int *change = INTEGER(CONDTARGETS);
		pModel->conditional(true);
		pModel->conditionalDependentVariable(CHAR(STRING_ELT(CONDVAR,0)));

		int i = 0;

		for (int group = 0; group < nGroups; group++)
		{
			Data * pData = (*pGroupData)[group];

			for (int period = 0;
				 period < pData->observationCount() - 1;
				 period++)
			{
				pModel->targetChange(pData, period, change[i]);
				i++;
			}
		}
	}
	/* get names vector for max degree */
	if (!Rf_isNull(MAXDEGREE))
	{
		SEXP Names = Rf_getAttrib(MAXDEGREE, R_NamesSymbol);

		for (int group = 0; group < nGroups; group++)
		{
			for (int i = 0; i < Rf_length(Names); i++)
			{
				Data * pData = (*pGroupData)[group];
				NetworkLongitudinalData * pNetworkData =
					pData->pNetworkData(CHAR(STRING_ELT(Names, i)));
				pNetworkData->maxDegree(INTEGER(MAXDEGREE)[i]);
			}
		}
	}
	/* get names vector for UniversalOffset */
	if (!Rf_isNull(UNIVERSALOFFSET))
	{
		SEXP Names = Rf_getAttrib(UNIVERSALOFFSET, R_NamesSymbol);

		for (int group = 0; group < nGroups; group++)
		{
			for (int i = 0; i < Rf_length(Names); i++)
			{
				Data * pData = (*pGroupData)[group];
				NetworkLongitudinalData * pNetworkData =
					pData->pNetworkData(CHAR(STRING_ELT(Names, i)));
				pNetworkData->universalOffset(REAL(UNIVERSALOFFSET)[i]);
			}
		}
	}
	/* set the parallel run flag on the model */
	if (!Rf_isNull(PARALLELRUN))
	{
		pModel->parallelRun(true);
	}
	/* get names vector for modeltype */
	if (!Rf_isNull(MODELTYPE))
	{
		SEXP Names = Rf_getAttrib(MODELTYPE, R_NamesSymbol);

		for (int group = 0; group < nGroups; group++)
		{
			for (int i = 0; i < Rf_length(Names); i++)
			{
				Data * pData = (*pGroupData)[group];
				NetworkLongitudinalData * pNetworkData =
					pData->pNetworkData(CHAR(STRING_ELT(Names, i)));
				pNetworkData->modelType(INTEGER(MODELTYPE)[i]);
			}
		}
	}
	/* get names vector for modeltype */
	if (!Rf_isNull(BEHMODELTYPE))
	{
		SEXP Names = Rf_getAttrib(BEHMODELTYPE, R_NamesSymbol);

		for (int group = 0; group < nGroups; group++)
		{
			for (int i = 0; i < Rf_length(Names); i++)
			{
				Data * pData = (*pGroupData)[group];
				BehaviorLongitudinalData * pBehaviorData =
					pData->pBehaviorData(CHAR(STRING_ELT(Names, i)));
				pBehaviorData->behModelType(INTEGER(BEHMODELTYPE)[i]);
			}
		}
	}
//	if (!Rf_isNull(MODELTYPE))
//	{
//		pModel->modelType(Rf_asInteger(MODELTYPE));
//	}
	// print out Data for profiling
	if (Rf_asInteger(PROFILEDATA))
	{
		printOutData((*pGroupData)[0]);
	}

	pModel->simpleRates(Rf_asInteger(SIMPLERATES));

	return R_NilValue;

}

SEXP getTargetActorStatistics(SEXP dataptr, SEXP modelptr, SEXP effectslist, SEXP parallelrun)
{
	vector<Data *> * pGroupData = (vector<Data *> *) R_ExternalPtrAddr(dataptr);
	Model * pModel = (Model *) R_ExternalPtrAddr(modelptr);

	if (!Rf_isNull(parallelrun))
	{
		pModel->parallelRun(true);
	}
	size_t nGroups = pGroupData->size();

	SEXP altStats = PROTECT(Rf_allocVector(VECSXP, nGroups));
	SEXP NETWORKTYPES = PROTECT(createRObjectAttributes(effectslist, altStats));
	int objEffects = Rf_length(NETWORKTYPES);

	for (size_t group = 0; group < nGroups; group++)
	{
		SET_VECTOR_ELT(altStats, group, Rf_allocVector(VECSXP, (*pGroupData)[group]->observationCount()));
		for (int p = 0; p < (*pGroupData)[group]->observationCount(); p++)
	{
			SET_VECTOR_ELT(VECTOR_ELT(altStats,group), p, Rf_allocVector(VECSXP, objEffects));
	}
	}

	for (size_t group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];
		for (int period = 0; period < pData->observationCount() - 1; period++)
	{
			State state(pData, period + 1);
			StatisticCalculator calculator(pData, pModel, &state, period, true);
			int actors = pData->rDependentVariableData()[0]->n();
			vector<double *> actorStatistics;
			getActorStatistics(effectslist, &calculator, &actorStatistics);
			for (unsigned e = 0; e < actorStatistics.size(); e++)
		 {
				SEXP actorStatsValues;
				PROTECT(actorStatsValues = Rf_allocVector(REALSXP,actors));
				double * astats = REAL(actorStatsValues);
				for (int i = 0; i < actors; i++)
			{
					astats[i] = actorStatistics.at(e)[i];
				}
				SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(altStats,group), period+1), e, actorStatsValues);
				UNPROTECT(1);
			}
		}
	}
	UNPROTECT(2);
	return altStats;
}

SEXP getTargetsChangeContributions(SEXP DATAPTR, SEXP MODELPTR, SEXP EFFECTSLIST, SEXP PARALLELRUN)
{
	vector<Data *> * pGroupData = (vector<Data *> *) R_ExternalPtrAddr(DATAPTR);
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	if (!Rf_isNull(PARALLELRUN))
	{
		pModel->parallelRun(true);
	}
	size_t nGroups = pGroupData->size();

	SEXP altStats = PROTECT(Rf_allocVector(VECSXP, nGroups));
	SEXP NETWORKTYPES = PROTECT(createRObjectAttributes(EFFECTSLIST, altStats));
	int objEffects = Rf_length(NETWORKTYPES);

	for (size_t group = 0; group < nGroups; group++)
	{
		SET_VECTOR_ELT(altStats, group, Rf_allocVector(VECSXP, (*pGroupData)[group]->observationCount()));
			for (int p = 0; p < (*pGroupData)[group]->observationCount(); p++)
		{
			SET_VECTOR_ELT(VECTOR_ELT(altStats,group), p, Rf_allocVector(VECSXP,objEffects));
		}
	}

	for (size_t group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];

		for (int period = 0; period < pData->observationCount() - 1; period++)
		{
			State State (pData, period + 1);
			// TODO altStats[0] == altStats[1] ??
			StatisticCalculator calculator(pData, pModel, &State, period, false, true);
			vector<vector<double * > > changeContributions;
			getChangeContributionStatistics(EFFECTSLIST, &calculator, &changeContributions);
			int actors = pData->rDependentVariableData()[0]->n();
			for(unsigned e = 0; e < changeContributions.size(); e++)
			{
				SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(altStats,group), period+1),
								e, Rf_allocVector(VECSXP,actors));
				int choices;
				if (strcmp(CHAR(STRING_ELT(NETWORKTYPES,e)), "behavior") == 0)
				{
					choices = 3;
				}
				else
				{
					choices = actors; // will not work for bipartite
					// But I did not find a convenient way to get knowledge of m() to here
				}
				for(int actor = 0; actor < actors; actor++)
				{
					SEXP actorsVal = PROTECT(Rf_allocVector(REALSXP, choices));
					double * d = REAL(actorsVal);
					for(int i = 0; i< Rf_length(actorsVal); i++)
					{
						d[i]=changeContributions.at(e).at(actor)[i];
					}
					SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(altStats, group), period+1), e), 
								actor, actorsVal);
					UNPROTECT(1);
				}
			}
		}
		State State (pData, 0);
		StatisticCalculator calculator(pData, pModel, &State, 0, false, true);
			vector<vector<double * > > changeContributions;
		getChangeContributionStatistics(EFFECTSLIST, &calculator, &changeContributions);
		int actors = pData->rDependentVariableData()[0]->n();
		for(unsigned e = 0; e < changeContributions.size(); e++)
		{
			SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(altStats, group), 0),
				e, Rf_allocVector(VECSXP, actors));
			int choices;
			if (strcmp(CHAR(STRING_ELT(NETWORKTYPES,e)), "behavior") == 0)
			{
				choices = 3;
			}
			else
			{						
					choices = actors; // will not work for bipartite
			}
			for(int actor = 0; actor < actors; actor++)
			{
				SEXP actorsVal = PROTECT(Rf_allocVector(REALSXP, choices));
				double * d = REAL(actorsVal);
				for(int i = 0; i< Rf_length(actorsVal); i++)
				{
					d[i]=changeContributions.at(e).at(actor)[i];
				}
				SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(altStats, group), 0), e), 
						actor, actorsVal);
				UNPROTECT(1);
			}
		}
	}
	UNPROTECT(2);
	return altStats;
}

/**
 *  Gets target values relative to the input data
 */
SEXP getTargets(SEXP DATAPTR, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP PARALLELRUN, SEXP RETURNACTORSTATISTICS,
	SEXP RETURNSTATICCHANGECONTRIBUTIONS)
{
	// select alternative result type
	int returnActorStatistics = sexp_to_int(RETURNACTORSTATISTICS, 0);
	int returnStaticChangeContributions = sexp_to_int(RETURNSTATICCHANGECONTRIBUTIONS, 0);
	if(returnActorStatistics + returnStaticChangeContributions >= 2)
	{
		Rf_error("returnActorStatistics and returnStaticChangeContributions are mutually exclusive");
	}
	if (returnActorStatistics) {
		return getTargetActorStatistics(DATAPTR, MODELPTR, EFFECTSLIST, PARALLELRUN);
	}
	if (returnStaticChangeContributions) {
		return getTargetsChangeContributions(DATAPTR, MODELPTR, EFFECTSLIST, PARALLELRUN);
	}

	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *) R_ExternalPtrAddr(DATAPTR);

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	if (!Rf_isNull(PARALLELRUN))
	{
		//TODO is this correct?
		pModel->parallelRun(true);
		// Or shouldn't it be:
		// 	if(Rf_asInteger(PARALLELRUN)==1)
		//  {
		//		pModel->parallelRun(true);
		//  }
	}
	size_t nGroups = pGroupData->size();
	int totObservations = totalPeriods(*pGroupData);

	// find the number of effects over all dependent variables:
	// sum of lengths of first columns:
	// for dimension of return vector
	int nEffects = 0;
	for (int i = 0; i < Rf_length(EFFECTSLIST); i++)
	{
		nEffects += Rf_length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}

	/* fra will contain the simulated statistics and must be initialised
	   to 0. Use rfra to reduce function evaluations. */
	SEXP fra = PROTECT(Rf_allocMatrix(REALSXP, nEffects, totObservations));
	double * rfra = REAL(fra);
	for (int i = 0; i < Rf_length(fra); i++)
	{
		rfra[i] = 0;
	}

	// find the targets: for each data object separately:
	// add them up on return to R (easier to check!)
	int periodFromStart = 0;

	for (size_t group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];

		for (int period = 0; period < pData->observationCount() - 1; period++)
		{
			periodFromStart++;
			//	EpochSimulation  Simulation(pData, pModel);
			//Simulation.initialize(period + 1);
			State State (pData, period + 1);
			//State State (&Simulation);
			StatisticCalculator Calculator (pData, pModel, &State,
				period);
			vector<double> statistic(nEffects);
			vector<double> score(nEffects); /* not used */
			vector<double> deriv(nEffects*nEffects); /* ABC not used */

			getStatistics(EFFECTSLIST, &Calculator, period,
					group, pData, (EpochSimulation *) 0, &statistic, &score); 
			//getStatistics(EFFECTSLIST, &Calculator, period,
			//			group, pData, &Simulation,
			//&statistic, &score);
			//	Rprintf("%f %f \n",statistic[1], statistic[2]);
			/* fill up matrices for  return value list */
			int iii = (periodFromStart - 1) * nEffects;
			for (unsigned effectNo = 0; effectNo < statistic.size(); effectNo++)
	{
				rfra[iii + effectNo] = statistic[effectNo];
			}
		}
	}
	UNPROTECT(1);
		return fra;
}

/**
 * Sets up a minimal chain and does pre burnin and burnin.
 * Processes a complete set of data objects, creating a chain for each
 * period and returning the address.
 */
/**
  * NOTE; FOR SOME CONFIGURATIONS OF STRUCTURAL ZEROS
 *  THIS RUNS INTO A HANG
 */
SEXP mlMakeChains(SEXP DATAPTR, SEXP MODELPTR,
		SEXP PROBS, SEXP PRMIN, SEXP PRMIB, SEXP MINIMUMPERM,
		SEXP MAXIMUMPERM, SEXP INITIALPERM, SEXP LOCALML)
{
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *) R_ExternalPtrAddr(DATAPTR);
	int nGroups = pGroupData->size();

	/* find total number of periods to process */
	int totObservations = totalPeriods(*pGroupData);

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);
	// create chain storage
	pModel->setupChainStore(totObservations);

	/* copy permutation lengths to the model */

	pModel->maximumPermutationLength(REAL(MAXIMUMPERM)[0]);
	pModel->minimumPermutationLength(REAL(MINIMUMPERM)[0]);
	pModel->initialPermutationLength(REAL(INITIALPERM)[0]);
	pModel->initializeCurrentPermutationLength();
	/* set MH probability values */
	pModel->insertDiagonalProbability(REAL(PROBS)[0]);
	pModel->cancelDiagonalProbability(REAL(PROBS)[1]);
	pModel->permuteProbability(REAL(PROBS)[2]);
	pModel->insertPermuteProbability(REAL(PROBS)[3]);
	pModel->deletePermuteProbability(REAL(PROBS)[4]);
	pModel->insertRandomMissingProbability(REAL(PROBS)[5]);
	pModel->deleteRandomMissingProbability(REAL(PROBS)[6]);
	//Rf_PrintValue(PROBS);

	double * prmin = REAL(PRMIN);
	double * prmib = REAL(PRMIB);

	SEXP minimalChains = PROTECT(Rf_allocVector(VECSXP, totObservations));
	SEXP currentChains = PROTECT(Rf_allocVector(VECSXP, totObservations));
	SEXP accepts = PROTECT(Rf_allocVector(VECSXP, totObservations));
	SEXP rejects = PROTECT(Rf_allocVector(VECSXP, totObservations));
	SEXP aborts = PROTECT(Rf_allocVector(VECSXP, totObservations));
	GetRNGstate();

	/* localML */
	int localML = 0;
	if (!Rf_isNull(LOCALML))
	{
		localML = Rf_asInteger(LOCALML);
	}
	pModel->localML(localML);

	int periodFromStart = 0;

	for (int group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];
		int observations = pData->observationCount() - 1;

		/* create the ML simulation object */
		MLSimulation * pMLSimulation = new MLSimulation(pData, pModel);

		pMLSimulation->simpleRates(pModel->simpleRates());

		for (int period = 0; period < observations; period ++)
		{
			// store for later on model
			pModel->missingNetworkProbability(prmin[periodFromStart]);
			pModel->missingBehaviorProbability(prmib[periodFromStart]);

			// put ones for this period on simulation object
			pMLSimulation->missingNetworkProbability(prmin[periodFromStart]);
			pMLSimulation->missingBehaviorProbability(prmib[periodFromStart]);

			pMLSimulation->currentPermutationLength(
				pModel->currentPermutationLength(period));

			/* initialize the chain: this also initializes the data */
			// does not initialise with previous period missing values yet
			pMLSimulation->pChain()->clear();
			pMLSimulation->connect(period);
			SEXP ch = PROTECT(getChainDFPlus(*(pMLSimulation->pChain()), true));
			SET_VECTOR_ELT(minimalChains, periodFromStart, ch);
			UNPROTECT(1);

			/* get the chain up to a reasonable length */
			pMLSimulation->preburnin();

			/* do some more steps */
			pMLSimulation->setUpProbabilityArray();

			int numSteps = 500;
			for (int i = 0; i < numSteps; i++)
			{
				pMLSimulation->MLStep();
			}

			/* store chain on Model after creating difference vectors */
			Chain * pChain = pMLSimulation->pChain();
			pMLSimulation->updateProbabilities(pChain,
				pChain->pFirst()->pNext(),
				pChain->pLast()->pPrevious());
			pChain->createInitialStateDifferences();
			pMLSimulation->createEndStateDifferences();
			pModel->chainStore(*pChain, periodFromStart);

			/* return chain as a list */
			SEXP ch1 = PROTECT(getChainList(*pChain));
			//PROTECT(ch1 = getChainDFPlus(*pChain, true));
			SET_VECTOR_ELT(currentChains, periodFromStart, ch1);
			UNPROTECT(1);

			/* get hold of the statistics for accept and reject */
			const vector < DependentVariable * > & rVariables =
				pMLSimulation->rVariables();
			int numberVariables = rVariables.size();

			SEXP accepts1 = PROTECT(Rf_allocMatrix(INTSXP, numberVariables, NBRTYPES));
			SEXP rejects1 = PROTECT(Rf_allocMatrix(INTSXP, numberVariables, NBRTYPES));
			SEXP aborts1 = PROTECT(Rf_allocVector(INTSXP, NBRTYPES));
			int * iaccepts = INTEGER(accepts1);
			int * irejects = INTEGER(rejects1);
			int * iaborts = INTEGER(aborts1);
			for (int i = 0; i < NBRTYPES; i++)
			{
				iaborts[i] = pMLSimulation->aborts(i);
				for (int j = 0; j < numberVariables; j++)
				{
					iaccepts[i + NBRTYPES * j] = rVariables[j]->acceptances(i);
					irejects[i + NBRTYPES * j] = rVariables[j]->rejections(i);
				}
			}
			SET_VECTOR_ELT(accepts, periodFromStart, accepts1);
			SET_VECTOR_ELT(rejects, periodFromStart, rejects1);
			SET_VECTOR_ELT(aborts, periodFromStart, aborts1);
			UNPROTECT(3);

			periodFromStart++;
			pModel->currentPermutationLength(period,
				pMLSimulation->currentPermutationLength());
		}
		delete pMLSimulation;
	}

	SEXP ans = PROTECT(Rf_allocVector(VECSXP, 5));
	SET_VECTOR_ELT(ans, 0, minimalChains);
	SET_VECTOR_ELT(ans, 1, currentChains);
	SET_VECTOR_ELT(ans, 2, accepts);
	SET_VECTOR_ELT(ans, 3, rejects);
	SET_VECTOR_ELT(ans, 4, aborts);

	PutRNGstate();

	UNPROTECT(6);
	return ans;
}

/**
 * Sets up chains in sub processes by copying them from input
 */
SEXP mlInitializeSubProcesses(SEXP DATAPTR, SEXP MODELPTR, SEXP PROBS,
		SEXP PRMIN, SEXP PRMIB, SEXP MINIMUMPERM, SEXP MAXIMUMPERM,
		SEXP INITIALPERM, SEXP CHAINS, SEXP LOCALML)
{
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *) R_ExternalPtrAddr(DATAPTR);

	int nGroups = pGroupData->size();

	/* find total number of periods to process */
	int totObservations = totalPeriods(*pGroupData);

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	// create chain storage
	pModel->setupChainStore(totObservations);

	/* copy permutation lengths to the model */

	pModel->maximumPermutationLength(REAL(MAXIMUMPERM)[0]);
	pModel->minimumPermutationLength(REAL(MINIMUMPERM)[0]);
	pModel->initialPermutationLength(REAL(INITIALPERM)[0]);
	pModel->initializeCurrentPermutationLength();
	/* set probability flags */
	pModel->insertDiagonalProbability(REAL(PROBS)[0]);
	pModel->cancelDiagonalProbability(REAL(PROBS)[1]);
	pModel->permuteProbability(REAL(PROBS)[2]);
	pModel->insertPermuteProbability(REAL(PROBS)[3]);
	pModel->deletePermuteProbability(REAL(PROBS)[4]);
	pModel->insertRandomMissingProbability(REAL(PROBS)[5]);
	pModel->deleteRandomMissingProbability(REAL(PROBS)[6]);

	double * prmin = REAL(PRMIN);
	double * prmib = REAL(PRMIB);

	int periodFromStart = 0;

	/* localML */
	int localML = sexp_to_int(LOCALML, 0);
	pModel->localML(localML);

	for (int group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];
		int observations = pData->observationCount() - 1;

		for (int period = 0; period < observations; period ++)
		{
			// store for later on model
			pModel->missingNetworkProbability(prmin[periodFromStart]);
			pModel->missingBehaviorProbability(prmib[periodFromStart]);

			/* copy the chain for this period onto the model */
			Chain * pChain = makeChainFromList(pData,
				VECTOR_ELT(CHAINS, periodFromStart), period);
			pModel->chainStore(*pChain, periodFromStart);

			periodFromStart++;
		}
	}

	return R_NilValue;
}
}

