/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07internals.cpp
 *
 * Description: This module contains routines used in to set up data in C.
 * Internal functions to C++: not visible from R.
 *****************************************************************************/
/**
 * @file
 * Internal routines used to set up the Data object with data from R.
 * Not visible from R.
 */
#include <vector>
#include <cstring>
// #include <Rinternals.h> // included by siena07internals.h
#include "siena07internals.h"
#include "data/Data.h"
#include "data/LongitudinalData.h"
#include "data/NetworkLongitudinalData.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ContinuousLongitudinalData.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/ConstantCovariate.h"
#include "data/ExogenousEvent.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "model/State.h"
#include "model/StatisticCalculator.h"
#include "data/ActorSet.h"
#include "model/EpochSimulation.h"
#include "model/SdeSimulation.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/MLSimulation.h"
#include "network/layers/PrimaryLayer.h"
#include "network/TieIterator.h"
#include "utils/Utils.h"
#include <Rinternals.h>


using namespace std;
using namespace siena;


/**
 * Matches column names with indices. The object Names is the names of the
 * effects data frame.
 */
void getColNos(SEXP Names, int * netTypeCol, int * nameCol, int * effectCol,
	int * parmCol, int * int1Col, int * int2Col, int * initValCol,
	int * typeCol, int * groupCol, int * periodCol, int * pointerCol,
	int * rateTypeCol, int * intptr1Col, int * intptr2Col, int * intptr3Col,
	int * settingCol)
{
	*netTypeCol = -1; /* net type */
	*nameCol = -1; /* network name */
	*effectCol = -1;  /* short name of effect */
	*parmCol = -1;
	*int1Col = -1;
	*int2Col = -1;
	*initValCol = -1;
	*typeCol = -1;
	*groupCol = -1;
	*periodCol = -1;
	*pointerCol = -1;
	*rateTypeCol = -1;
	*intptr1Col = -1;
	*intptr2Col = -1;
	*intptr3Col = -1;
	*settingCol = -1;

	int n = Rf_length(Names);
	for (int j = 0; j < n; j++)
	{
		if (strcmp(CHAR(STRING_ELT(Names, j)), "netType") == 0)
		{
			*netTypeCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "name") == 0)
		{
			*nameCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "shortName") == 0)
		{
			*effectCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "parm") == 0)
		{
			*parmCol = j;
		}

		if (strcmp(CHAR(STRING_ELT(Names, j)), "interaction1") == 0)
		{
			*int1Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "interaction2") == 0)
		{
			*int2Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "initialValue") == 0)
		{
			*initValCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "type") == 0)
		{
			*typeCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "group") == 0)
		{
			*groupCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "period") == 0)
		{
			*periodCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "effectPtr") == 0)
		{
			*pointerCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "rateType") == 0)
		{
			*rateTypeCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "effect1") == 0)
		{
			*intptr1Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "effect2") == 0)
		{
			*intptr2Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "effect3") == 0)
		{
			*intptr3Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "setting") == 0)
		{
			*settingCol = j;
		}
	}
	if (*netTypeCol < 0)
	{
		Rf_error("cannot find nettype");
	}
	if (*nameCol < 0)
	{
		Rf_error("cannot find network name");
	}
	if (*effectCol < 0)
	{
		Rf_error("cannot find effectName");
	}
	if (*parmCol < 0)
	{
		Rf_error("cannot find internal parameter");
	}
	if (*int1Col < 0)
	{
		Rf_error("cannot find interaction1");
        }

	if (*int2Col < 0)
	{
		Rf_error("cannot find interaction2");
	}
	if (*initValCol < 0)
	{
		Rf_error("cannot find initial value");
	}
	if (*groupCol < 0)
	{
		Rf_error("cannot find group");
	}
	if (*periodCol < 0)
	{
		Rf_error("cannot find period");
	}
	if (*pointerCol < 0)
	{
		Rf_error("cannot find effect pointer");
	}
	if (*rateTypeCol < 0)
	{
		Rf_error("cannot find rate type");
	}
	if (*intptr1Col < 0)
	{
		Rf_error("cannot find effect1");
	}
	if (*intptr2Col < 0)
	{
		Rf_error("cannot find effect2");
	}
	if (*intptr3Col < 0)
	{
		Rf_error("cannot find effect3");
	}
	if (*settingCol < 0)
	{
		Rf_error("cannot find setting col; reconstruct effects object with this version of RSiena");
	}
	//Rprintf("%d parmcol\n", *parmCol);
}

/**
 *  updates the parameter values for each of the effects.
 */
void updateParameters(SEXP EFFECTSLIST, SEXP THETA, vector<Data *> *
		pGroupData, Model * pModel)
{
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

	int thetasub = -1;
	/* find each effect and update its weight */
	for (int net = 0; net < Rf_length(EFFECTSLIST); net++)
	{
		const char * networkName =
			CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, net),
									   nameCol), 0));
		SEXP EFFECTS = VECTOR_ELT(EFFECTSLIST, net);
		for (int eff = 0; eff < Rf_length(VECTOR_ELT(EFFECTS,0)); eff++)
		{
			thetasub = thetasub + 1;
			const char * effectName =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol),  eff));
			double currentValue = REAL(THETA)[thetasub];
			const char * effectType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), eff));
			const char * netType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, netTypeCol), eff));
			const char * setting =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, settingCol), eff));

			if (strcmp(effectType, "rate") == 0 &&
				strcmp(effectName, "Rate") == 0)
			{
				int group =
					INTEGER(VECTOR_ELT(EFFECTS, groupCol))[eff]  - 1;
				int period =
					INTEGER(VECTOR_ELT(EFFECTS, periodCol))[eff] - 1;
				/* find the network */
				Data * pData = (*pGroupData)[group];

				if (strcmp(setting, "") == 0)
				{
					if (strcmp(netType, "behavior") == 0)
					{
						BehaviorLongitudinalData * pNetwork =
							pData->pBehaviorData(networkName);
						pModel->basicRateParameter(pNetwork,
							period,
							currentValue);
					}
					else
					{
						NetworkLongitudinalData * pNetwork =
							pData->pNetworkData(networkName);
						pModel->basicRateParameter(pNetwork,
							period, currentValue);
					}
				}
				else
				{
				 	if (!(strcmp(netType, "behavior") == 0))
				 	{
				 		NetworkLongitudinalData * pNetwork =
				 			pData->pNetworkData(networkName);
				 		pModel->settingRateParameter(pNetwork, setting, period,
				 			currentValue);
				 	}
				 	else
				 	{
				 		Rf_error("setting found for behavior variable %s",
				 			networkName);
				 	}
				}
			}
			else if (strcmp(effectType, "rate") == 0 &&
					 strcmp(effectName, "scale") == 0)
			{
				int period = INTEGER(VECTOR_ELT(EFFECTS, periodCol))[eff] - 1;
				if (strcmp(setting, "") == 0)
				{
					pModel->basicScaleParameter(period, currentValue);
				}
			else
			{
					Rf_error("setting found for behavior variable %s", 
						networkName);
				}
			}
			else // no rate or scale effect
			{
				EffectInfo * pEffectInfo =
					(EffectInfo *) R_ExternalPtrAddr(
						VECTOR_ELT(VECTOR_ELT(EFFECTS, pointerCol), eff));
				pEffectInfo->parameter(currentValue);
			}
		}
	}

	UNPROTECT(1);
	return;
}


/**
 * Create one observation for a one mode Network: ties, missing, structural
 *
 */
void setupOneModeNetwork(SEXP ONEMODE,
	OneModeNetworkLongitudinalData * pNetworkData,
	int observation)
{
	/* one mode networks are passed in as list of edgelists with attributes
	 giving the size of the network - not checked yet*/

	// Tie values
	//Rprintf("%x\n", pNetworkData);
	SEXP ONEMODEVALS = VECTOR_ELT(ONEMODE, 0);
	int *start = INTEGER(ONEMODEVALS);
	int listlen = Rf_ncols(ONEMODEVALS);
	int pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->tieValue(i - 1, j - 1, observation, val);
	}

	// Missingness

	ONEMODEVALS = VECTOR_ELT(ONEMODE, 1);
	start = INTEGER(ONEMODEVALS);
	listlen = Rf_ncols(ONEMODEVALS);
	pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->missing(i - 1, j - 1, observation, val);
	}

	// Structural ties

	ONEMODEVALS = VECTOR_ELT(ONEMODE, 2);
	start = INTEGER(ONEMODEVALS);
	listlen = Rf_ncols(ONEMODEVALS);
	pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->structural(i - 1, j - 1, observation, val);
	}
}


/**
 * Create all observations for a one mode Network
 *
 */
void setupOneModeObservations(const std::string& name, SEXP ONEMODES,
		OneModeNetworkLongitudinalData * pOneModeNetworkLongitudinalData)
{
	int observations = Rf_length(ONEMODES);
	if (observations != pOneModeNetworkLongitudinalData->observationCount())
	{
		Rf_error("wrong number of observations for one-mode network");
	}
	SEXP uo;
	PROTECT(uo = Rf_install("uponly"));
	SEXP uponly = Rf_getAttrib(ONEMODES, uo);
	SEXP dow;
	PROTECT(dow = Rf_install("downonly"));
	SEXP downonly = Rf_getAttrib(ONEMODES, dow);

	for (int period = 0; period < (observations - 1); period++)
	{
		pOneModeNetworkLongitudinalData->upOnly(period,
				LOGICAL(uponly)[period]);
		pOneModeNetworkLongitudinalData->downOnly(period,
				LOGICAL(downonly)[period]);
	}
	for (int period = 0; period < observations; period++)
	{
		setupOneModeNetwork(VECTOR_ELT(ONEMODES, period),
				pOneModeNetworkLongitudinalData, period);
	}
	UNPROTECT(2);
}

/**
 * Create one group of one mode Networks
 *
 */
void setupOneModeGroup(SEXP ONEMODEGROUP, Data * pData)
{
	int nOneMode = Rf_length(ONEMODEGROUP);

	for (int oneMode = 0; oneMode < nOneMode; oneMode++)
	{
		SEXP ONEMODES = VECTOR_ELT(ONEMODEGROUP, oneMode);
		SEXP as = PROTECT(Rf_install("nodeSet"));
		SEXP actorSet = PROTECT(Rf_getAttrib(ONEMODES, as));
		SEXP symm = PROTECT(Rf_install("symmetric"));
		SEXP symmetric = PROTECT(Rf_getAttrib(ONEMODES, symm));
		SEXP balm = PROTECT(Rf_install("balmean"));
		SEXP balmean = PROTECT(Rf_getAttrib(ONEMODES, balm));
		SEXP strm = PROTECT(Rf_install("structmean"));
		SEXP structmean = PROTECT(Rf_getAttrib(ONEMODES, strm));
		SEXP avin = PROTECT(Rf_install("averageInDegree"));
		SEXP averageInDegree = PROTECT(Rf_getAttrib(ONEMODES, avin));
		SEXP avout = PROTECT(Rf_install("averageOutDegree"));
		SEXP averageOutDegree = PROTECT(Rf_getAttrib(ONEMODES, avout));

		SEXP nm = PROTECT(Rf_install("name"));
		SEXP name = Rf_getAttrib(ONEMODES, nm);
		const ActorSet* pActorSet = pData->pActorSet(CHAR(STRING_ELT(actorSet, 0)));
		const char* cname = CHAR(STRING_ELT(name, 0));
		OneModeNetworkLongitudinalData *  pNetData = pData->createOneModeNetworkData(cname, pActorSet);

		// parse settings
		SEXP settingsSymbol = PROTECT(Rf_install("settingsinfo"));
		SEXP settingsList = PROTECT(Rf_getAttrib(ONEMODES, settingsSymbol));
		for (int j = 0; j < Rf_length(settingsList); j++)
		{
			SEXP settingInfo = VECTOR_ELT(settingsList, j);
			SEXP infoNames = Rf_getAttrib(settingInfo, R_NamesSymbol);
			std::string id, type, covar, only;
//			Rprintf("setting %d\n", j);
			// parse key value list
			for (int k = 0; k < Rf_length(settingInfo); k++) {
				// Rprintf("key value %d\n", k);
				const char* key = CHAR(STRING_ELT(infoNames, k));
				const char* value = CHAR(STRING_ELT(VECTOR_ELT(settingInfo, k), 0));
				// Rprintf("%s = %s\n", key, value);
				if      (strcmp(key, "id") == 0)    id    = std::string(value);
				else if (strcmp(key, "type") == 0)  type  = std::string(value);
				else if (strcmp(key, "covar") == 0) covar = std::string(value);
				else if (strcmp(key, "only") == 0)  only  = std::string(value);
			}
			// parse strings to C++ types
			Permission_Type permType = Permission_Type::BOTH;
			if      (only == "up")   permType = Permission_Type::UP;
			else if (only == "down") permType = Permission_Type::DOWN;
			// asserts
			if (id.length() == 0) Rf_error("settings id should not be empty");
			if (type.length() == 0) Rf_error("settings type should not be empty");
			// add it
			Rprintf("%s %s %s %s\n", id.c_str(), type.c_str(), covar.c_str(), only.c_str());
			pNetData->addSettingName(id, type, covar, permType);
		}

		// check first two setting types
		if (pNetData->rSettingNames().size() == 1) {
			Rf_error("if setting are present use universal and primary");
		} else if (pNetData->rSettingNames().size() >= 2) {
			if (pNetData->rSettingNames().at(0).getSettingType() != "universal")
				Rf_error("first setting should be type=universal");
			if (pNetData->rSettingNames().at(1).getSettingType() != "primary")
				Rf_error("second setting should be type=primary");
		}

		pNetData->symmetric(*(LOGICAL(symmetric)));
		pNetData->balanceMean(*(REAL(balmean)));
		pNetData->structuralMean(*(REAL(structmean)));
		pNetData->averageInDegree(*(REAL(averageInDegree)));
		pNetData->averageOutDegree(*(REAL(averageOutDegree)));
		setupOneModeObservations(cname, ONEMODES, pNetData);
		// Once all network data has been stored, calculate some statistical
		// properties of that data.
		pNetData->calculateProperties();

		// add the primary setting for each observation
		if (pNetData->rSettingNames().size() >= 1) {
			std::string settingName = "primary(" + std::string(cname) + ")";
			OneModeNetworkLongitudinalData* pSettingNetwork =
				pData->createOneModeSimNetworkData(settingName.c_str(), pActorSet);

			// Copy the primary layer to the data objects.  Cannot simply use
			// `setupOneModeObservations` (with an object from R) because that
			// expects a different format.
			for (int obs = 0; obs < pNetData->observationCount(); obs++) {
				const Network& net = *pNetData->pNetwork(obs);
				PrimaryLayer l;
				l.onInitializationEvent(net);
				for (TieIterator tie = l.pLayer()->ties(); tie.valid(); tie.next()) {
					pSettingNetwork->tieValue(tie.ego(), tie.alter(), obs, 1);
					// TODO ignore missing and structural?
					pSettingNetwork->missing(tie.ego(), tie.alter(), obs, 0);
					pSettingNetwork->structural(tie.ego(), tie.alter(), obs, 0);
				}
				l.onNetworkDisposeEvent(net);
			}
			pSettingNetwork->calculateProperties();
		}

		UNPROTECT(15);
	}
}

/**
 * Create one observation for a bipartite Network: ties, missing, structural
 *
 */
void setupBipartiteNetwork(SEXP BIPARTITE,
	NetworkLongitudinalData * pNetworkData,
	int observation)
{
	/* bipartite networks are passed in as list of edgelists with attributes
	 giving the size of the network - not checked yet*/

	// Tie values

	SEXP BIPARTITEVALS = VECTOR_ELT(BIPARTITE, 0);
	int *start = INTEGER(BIPARTITEVALS);
	int listlen = Rf_ncols(BIPARTITEVALS);
	int pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->tieValue(i - 1, j - 1, observation, val);
	}

	// Missingness

	BIPARTITEVALS = VECTOR_ELT(BIPARTITE, 1);
	start = INTEGER(BIPARTITEVALS);
	listlen = Rf_ncols(BIPARTITEVALS);
	pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->missing(i - 1, j - 1, observation, val);
	}

	// Structural ties

	BIPARTITEVALS = VECTOR_ELT(BIPARTITE, 2);
	start = INTEGER(BIPARTITEVALS);
	listlen = Rf_ncols(BIPARTITEVALS);
	pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->structural(i - 1, j - 1, observation, val);
	}
}


/**
 * Create all observations for a bipartite Network
 *
 */
void setupBipartiteObservations(SEXP BIPARTITES,
	NetworkLongitudinalData *
	pNetworkLongitudinalData)

{
    int observations = Rf_length(BIPARTITES);
    if (observations != pNetworkLongitudinalData->observationCount())
    {
		Rf_error ("wrong number of observations in bipartite");
    }
    SEXP uo;
    PROTECT(uo = Rf_install("uponly"));
    SEXP uponly = Rf_getAttrib(BIPARTITES, uo);
    SEXP dow;
    PROTECT(dow = Rf_install("downonly"));
    SEXP downonly = Rf_getAttrib(BIPARTITES, dow);
    for (int period = 0; period < (observations - 1); period++)
    {
        pNetworkLongitudinalData->upOnly(period,
			LOGICAL(uponly)[period]);
        pNetworkLongitudinalData->downOnly(period,
			LOGICAL(downonly)[period]);
    }
    for (int period = 0; period < observations; period++)
    {
    	setupBipartiteNetwork(VECTOR_ELT(BIPARTITES, period),
			pNetworkLongitudinalData,
			period);
    }
    UNPROTECT(2);
}
/**
 * Create one group of bipartite Networks
 *
 */
void setupBipartiteGroup(SEXP BIPARTITEGROUP, Data * pData)
{
	int nBipartite = Rf_length(BIPARTITEGROUP);

	for (int bipartite = 0; bipartite < nBipartite; bipartite++)
	{
		SEXP as;
		PROTECT(as = Rf_install("nodeSet"));
		SEXP actorSet = Rf_getAttrib(VECTOR_ELT(BIPARTITEGROUP, bipartite), as);
		SEXP nm;
		PROTECT(nm = Rf_install("name"));
		SEXP name = Rf_getAttrib(VECTOR_ELT(BIPARTITEGROUP, bipartite), nm);
		SEXP avout;
		PROTECT(avout = Rf_install("averageOutDegree"));
		SEXP averageOutDegree = Rf_getAttrib(VECTOR_ELT(BIPARTITEGROUP,
					bipartite), avout);
		const ActorSet * pSenders = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 0)));
		const ActorSet * pReceivers = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 1)));
		NetworkLongitudinalData *  pNetworkLongitudinalData =
			pData->createNetworkData(CHAR(STRING_ELT(name, 0)),
					pSenders, pReceivers);
		pNetworkLongitudinalData->averageOutDegree(*(REAL(averageOutDegree)));
		setupBipartiteObservations(VECTOR_ELT(BIPARTITEGROUP, bipartite),
				pNetworkLongitudinalData);

		// Once all network data has been stored, calculate some
		// statistical properties of that data.

		pNetworkLongitudinalData->calculateProperties();
		UNPROTECT(3);
	}
}

/**
 * Create all observations for a behavior Network
 *
 */
void setupBehavior(SEXP BEHAVIOR, BehaviorLongitudinalData * pBehaviorData)
{
	int observations = Rf_ncols(VECTOR_ELT(BEHAVIOR, 0));

	if (observations != pBehaviorData->observationCount())
	{
		Rf_error ("wrong number of observations in Behavior");
	}
	int nActors = Rf_nrows(VECTOR_ELT(BEHAVIOR, 0));

	if (nActors != pBehaviorData->n())
	{
		Rf_error ("wrong number of actors");
	}
	int * start = INTEGER(VECTOR_ELT(BEHAVIOR, 0));
	int * missing = LOGICAL(VECTOR_ELT(BEHAVIOR, 1));

	for (int period = 0; period < observations; period++)
	{
		for (int actor = 0; actor < nActors; actor++)
		{
			pBehaviorData->value(period, actor, *start++);
			pBehaviorData->missing(period, actor, *missing++);
		}
	}
	SEXP uo;
	PROTECT(uo = Rf_install("uponly"));
	SEXP uponly = Rf_getAttrib(VECTOR_ELT(BEHAVIOR, 0), uo);
	SEXP dow;
	PROTECT(dow = Rf_install("downonly"));
	SEXP downonly = Rf_getAttrib(VECTOR_ELT(BEHAVIOR,0), dow);
	for (int period = 0; period < (observations - 1); period++)
	{
		pBehaviorData->upOnly(period, LOGICAL(uponly)[period]);
		pBehaviorData->downOnly(period, LOGICAL(downonly)[period]);
	}
	SEXP sim;
	PROTECT(sim = Rf_install("simMean"));
	SEXP simMean = Rf_getAttrib(VECTOR_ELT(BEHAVIOR,0), sim);
	pBehaviorData->similarityMean(REAL(simMean)[0]);
	SEXP sims;
	PROTECT(sims = Rf_install("simMeans"));
	SEXP simMeans = Rf_getAttrib(VECTOR_ELT(BEHAVIOR, 0), sims);
	SEXP simNames;
	PROTECT(simNames = Rf_getAttrib(simMeans, R_NamesSymbol));
	int numberNetworks = Rf_length(simMeans);
	for (int net = 0; net < numberNetworks; net++)
	{
		pBehaviorData->similarityMeans(REAL(simMeans)[net],
				CHAR(STRING_ELT(simNames, net)));
	}

	// Now that the values are set, calculate some important statistics
	pBehaviorData->calculateProperties();
	UNPROTECT(5);
}

/**
 * Create one group of Behavior Networks
 *
 */
void setupBehaviorGroup(SEXP BEHGROUP, Data *pData)
{
	int nBehavior = Rf_length(BEHGROUP);

	for (int behavior= 0; behavior < nBehavior; behavior++)
	{
		SEXP as;
		PROTECT(as = Rf_install("nodeSet"));
		SEXP actorSet = Rf_getAttrib(VECTOR_ELT(VECTOR_ELT(BEHGROUP, behavior), 0),
				as);

		SEXP nm;
		PROTECT(nm = Rf_install("name"));
		SEXP name = Rf_getAttrib(VECTOR_ELT(VECTOR_ELT(BEHGROUP, behavior), 0),
				nm);

		const ActorSet * pActorSet = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 0)));
		BehaviorLongitudinalData * pBehaviorData =
			pData->createBehaviorData(CHAR(STRING_ELT(name, 0)), pActorSet);
		//	Rprintf("%x\n", pBehaviorData);
		setupBehavior(VECTOR_ELT(BEHGROUP, behavior), pBehaviorData);
		UNPROTECT(2);
	}
}

/**
 * Create all observations for a continuous dependent variable
 *
 */
void setupContinuous(SEXP CONTINUOUS, ContinuousLongitudinalData *
	pContinuousData)
{
    int observations = Rf_ncols(VECTOR_ELT(CONTINUOUS, 0));

    if (observations != pContinuousData->observationCount())
    {
		Rf_error ("wrong number of observations in Continuous");
    }
    int nActors = Rf_nrows(VECTOR_ELT(CONTINUOUS, 0));

    if (nActors != pContinuousData->n())
    {
        Rf_error ("wrong number of actors");
    }
    double * start = REAL(VECTOR_ELT(CONTINUOUS, 0));
	int * missing = LOGICAL(VECTOR_ELT(CONTINUOUS, 1));

    for (int period = 0; period < observations; period++)
    {
        for (int actor = 0; actor < nActors; actor++)
        {
			pContinuousData->value(period, actor, *start++);
			pContinuousData->missing(period, actor, *missing++);
        }
    }
    SEXP uo;
    PROTECT(uo = Rf_install("uponly"));
    SEXP uponly = Rf_getAttrib(VECTOR_ELT(CONTINUOUS, 0), uo);
    SEXP dow;
    PROTECT(dow = Rf_install("downonly"));
    SEXP downonly = Rf_getAttrib(VECTOR_ELT(CONTINUOUS,0), dow);
    for (int period = 0; period < (observations - 1); period++)
    {
        pContinuousData->upOnly(period, LOGICAL(uponly)[period]);
        pContinuousData->downOnly(period, LOGICAL(downonly)[period]);
    }
    SEXP sim;
    PROTECT(sim = Rf_install("simMean"));
    SEXP simMean = Rf_getAttrib(VECTOR_ELT(CONTINUOUS,0), sim);
	pContinuousData->similarityMean(REAL(simMean)[0]);
	SEXP sims;
	PROTECT(sims = Rf_install("simMeans"));
	SEXP simMeans = Rf_getAttrib(VECTOR_ELT(CONTINUOUS, 0), sims);
	SEXP simNames;
	PROTECT(simNames = Rf_getAttrib(simMeans, R_NamesSymbol));
	int numberNetworks = Rf_length(simMeans);
	for (int net = 0; net < numberNetworks; net++)
	{
		pContinuousData->similarityMeans(REAL(simMeans)[net],
			CHAR(STRING_ELT(simNames, net)));
	}

    // Now that the values are set, calculate some important statistics
	pContinuousData->calculateProperties();
	UNPROTECT(5);
}
/**
 * Create one group of Continuous dependent variables
 *
 */
void setupContinuousGroup(SEXP CONTGROUP, Data *pData)
{
    int nCont = Rf_length(CONTGROUP);

    for (int continuous = 0; continuous < nCont; continuous++)
    {
		SEXP as;
		PROTECT(as = Rf_install("nodeSet"));
        SEXP actorSet = Rf_getAttrib(VECTOR_ELT(VECTOR_ELT(CONTGROUP, continuous), 0),
								  as);

        SEXP nm;
        PROTECT(nm = Rf_install("name"));
        SEXP name = Rf_getAttrib(VECTOR_ELT(VECTOR_ELT(CONTGROUP, continuous), 0),
							  nm);

        const ActorSet * myActorSet = pData->pActorSet(CHAR(STRING_ELT(
                                                                actorSet, 0)));
		ContinuousLongitudinalData * pContinuousData =
			pData->createContinuousData(CHAR(STRING_ELT(name, 0)), myActorSet);
		setupContinuous(VECTOR_ELT(CONTGROUP, continuous), pContinuousData);
        UNPROTECT(2);
    }
}
/**
 * Create a constant covariate
 *
 */
void setupConstantCovariate(SEXP COCOVAR,
		ConstantCovariate * pConstantCovariate)
{
	int nActors = Rf_length(COCOVAR);
	// Rprintf("%x\n", pConstantCovariate);
	double * start = REAL(COCOVAR);
	SEXP mn;
	PROTECT(mn = Rf_install("mean"));
	SEXP ans = Rf_getAttrib(COCOVAR, mn);
	double mean = REAL(ans)[0];
	SEXP cn;
	PROTECT(cn = Rf_install("centered"));
	ans = Rf_getAttrib(COCOVAR, cn);
	bool centered = LOGICAL(ans)[0];

	// extract imputationValues if provided by user
	SEXP im;
	PROTECT(im = Rf_install("imputationValues"));
	ans = Rf_getAttrib(COCOVAR, im);
	bool impute = FALSE;
	double * imputationValues;
	if(!Rf_isNull(ans))
	{
		impute = TRUE;
		imputationValues = REAL(ans);
		//		Rprintf("We have something to impute\n");
	}

	for (int actor = 0; actor < nActors; actor++)
	{
		double value = *start++;
		if (ISNAN(value))
		{
			if (impute) // imputationValues already were centered, if necessary
			{
				pConstantCovariate->value(actor, imputationValues[actor]);
				//				Rprintf("We impute value %f for actor %d\n",
				//					imputationValues[actor], actor + 1);
			}
			else if (centered) // no user input provided
			{
				pConstantCovariate->value(actor, 0);
				//				Rprintf("We use 0 for actor %d\n", actor + 1);
			}
			else // no user input provided, not centered
			{
				pConstantCovariate->value(actor, mean);
				//				Rprintf("We use the mean %f for actor %d\n", mean, actor + 1);
			}
			pConstantCovariate->missing(actor, 1);
		}
		else
		{
			pConstantCovariate->value(actor, value);
			pConstantCovariate->missing(actor, 0);

		}
	}
	UNPROTECT(3);
}
/**
 * Create one group of constant covariates
 *
 */
void setupConstantCovariateGroup(SEXP COCOVARGROUP, Data *pData)
{
	int nConstantCovariate = Rf_length(COCOVARGROUP);
	//    Rprintf("nConstantCovariate %d\n", nConstantCovariate);
	for (int constantCovariate = 0; constantCovariate < nConstantCovariate;
			constantCovariate++)
	{
		SEXP as;
		PROTECT(as = Rf_install("nodeSet"));
		SEXP actorSet = Rf_getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
				as);
		SEXP nm;
		PROTECT(nm = Rf_install("name"));
		SEXP name = Rf_getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate), nm);
		const ActorSet * pActorSet = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 0)));
		int nActors = Rf_length(VECTOR_ELT(COCOVARGROUP, constantCovariate));
		//    Rprintf("nactors %d\n", nActors);

		if (nActors != pActorSet->n())
		{
			Rf_error ("wrong number of actors");
		}
		ConstantCovariate * pConstantCovariate =
			pData->createConstantCovariate(CHAR(STRING_ELT(name, 0)),
					pActorSet);
		setupConstantCovariate(VECTOR_ELT(COCOVARGROUP,	constantCovariate),
				pConstantCovariate);
		SEXP mn;
		PROTECT(mn = Rf_install("mean"));
		SEXP obsmean = Rf_getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate), mn);
		SEXP cn;
		PROTECT(cn = Rf_install("centered"));
		SEXP ans = Rf_getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate), cn);
		bool centered = LOGICAL(ans)[0];
		if (centered)
		{
			pConstantCovariate->mean(0);
		}
		else
		{
			pConstantCovariate->mean(REAL(obsmean)[0]);
		}
		SEXP sim;
		PROTECT(sim = Rf_install("simMean"));
		SEXP simMean = Rf_getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
				sim);
		pConstantCovariate->similarityMean(REAL(simMean)[0]);
		SEXP sims;
		PROTECT(sims = Rf_install("simMeans"));
		SEXP simMeans = Rf_getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
				sims);
		SEXP simNames;
		PROTECT(simNames = Rf_getAttrib(simMeans, R_NamesSymbol));
		int numberNetworks = Rf_length(simMeans);
		for (int net = 0; net < numberNetworks; net++)
		{
			pConstantCovariate->similarityMeans(REAL(simMean)[net],
					CHAR(STRING_ELT(simNames, net)));
		}
		SEXP range;
		PROTECT(range = Rf_install("range"));
		SEXP Range = Rf_getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
				range);
		pConstantCovariate->range(REAL(Range)[0]);
		UNPROTECT(8);
	}
}

/**
 * Create all observations for a changing covariate
 *
 */
void setupChangingCovariate(SEXP VARCOVAR,
		ChangingCovariate * pChangingCovariate)
{
	int observations = Rf_ncols(VARCOVAR);
	int nActors = Rf_nrows(VARCOVAR);
	double * start = REAL(VARCOVAR);
	SEXP mn;
	PROTECT(mn = Rf_install("mean"));
	SEXP ans = Rf_getAttrib(VARCOVAR, mn);
	double mean = REAL(ans)[0];
	SEXP cn;
	PROTECT(cn = Rf_install("centered"));
	ans = Rf_getAttrib(VARCOVAR, cn);
	bool centered = LOGICAL(ans)[0];

	// extract imputationValues if provided by user
	SEXP im;
	PROTECT(im = Rf_install("imputationValues"));
	ans = Rf_getAttrib(VARCOVAR, im);
	bool impute = FALSE;
	double * imputationValues = 0;
	if(!Rf_isNull(ans))
	{
		impute = TRUE;
		imputationValues = REAL(ans);
	}

	for (int period = 0; period < observations; period++)
	{
		for (int actor = 0; actor < nActors; actor++)
		{
			double value = *start++;
			double imputationValue;
			if (impute)
			{
				imputationValue = *imputationValues++;
			}

			if (ISNAN(value))
			{
				if (impute) // imputationValues have been centered, if necessary
				{
					pChangingCovariate->value(actor, period, imputationValue);
					//					Rprintf("We impute value %f for actor %d in period %d\n",
					//										imputationValue, actor + 1, period + 1);
				}
				else if (centered) // no user input provided
				{
					pChangingCovariate->value(actor, period, 0);
					//					Rprintf("We use 0 for actor %d in period %d\n",
					//								actor + 1, period + 1);
				}
				else // no user input provided, not centered
				{
					pChangingCovariate->value(actor, period, mean);
					//					Rprintf("We use the mean %f for actor %d in period %d\n",
					//								mean, actor + 1, period + 1);
				}
				pChangingCovariate->missing(actor, period, 1);
			}
			else
			{
				pChangingCovariate->value(actor, period, value);
				pChangingCovariate->missing(actor, period, 0);
			}
		}
	}
	UNPROTECT(3);
}

/**
 * Create one group of changing covariates
 *
 */
void setupChangingCovariateGroup(SEXP VARCOVARGROUP, Data *pData)
{
	if (Rf_length(VARCOVARGROUP) == 0)
		return;
	int observations = Rf_ncols(VECTOR_ELT(VARCOVARGROUP,0));
	if (observations != pData->observationCount() - 1)
	{
		Rf_error ("wrong number of observations in Changing Covariate");
	}
	int nChangingCovariate = Rf_length(VARCOVARGROUP);
	for (int changingCovariate = 0;
			changingCovariate < nChangingCovariate;
			changingCovariate++)
	{
		SEXP as;
		PROTECT(as = Rf_install("nodeSet"));
		SEXP actorSet = Rf_getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
				as);
		SEXP nm;
		PROTECT(nm = Rf_install("name"));
		SEXP name = Rf_getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
				nm);
		const ActorSet * pActorSet = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 0)));
		int nActors = Rf_nrows(VECTOR_ELT(VARCOVARGROUP,changingCovariate));

		if (nActors != pActorSet->n())
		{
			Rf_error ("wrong number of actors");
		}
		ChangingCovariate * pChangingCovariate =
			pData->createChangingCovariate(CHAR(STRING_ELT(name, 0)),
					pActorSet);
		setupChangingCovariate(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
				pChangingCovariate);
		SEXP mn;
		PROTECT(mn = Rf_install("mean"));
		SEXP obsmean = Rf_getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate), mn);
		SEXP cn;
		PROTECT(cn = Rf_install("centered"));
		SEXP ans = Rf_getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate), cn);
		bool centered = LOGICAL(ans)[0];
		if (centered)
		{
			pChangingCovariate->mean(0);
		}
		else
		{
			pChangingCovariate->mean(REAL(obsmean)[0]);
		}
		SEXP sim;
		PROTECT(sim = Rf_install("simMean"));
		SEXP simMean = Rf_getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
				sim);
		pChangingCovariate->similarityMean(REAL(simMean)[0]);
		SEXP sims;
		PROTECT(sims = Rf_install("simMeans"));
		SEXP simMeans = Rf_getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
				sims);
		SEXP simNames;
		PROTECT(simNames = Rf_getAttrib(simMeans, R_NamesSymbol));
		int numberNetworks = Rf_length(simMeans);
		for (int net = 0; net < numberNetworks; net++)
		{
			pChangingCovariate->similarityMeans(REAL(simMean)[net],
					CHAR(STRING_ELT(simNames, net)));
		}
		SEXP range;
		PROTECT(range = Rf_install("range"));
		SEXP Range = Rf_getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
				range);
		pChangingCovariate->range(REAL(Range)[0]);
		UNPROTECT(8);
	}
}

/**
 * Create a constant dyadic covariate
 *
 */
void setupDyadicCovariate(SEXP DYADVAR,
		ConstantDyadicCovariate * pConstantDyadicCovariate)
{
	double *start = REAL(VECTOR_ELT(DYADVAR, 0));
	double *missingstart = REAL(VECTOR_ELT(DYADVAR, 1));
	int listlen = Rf_ncols(VECTOR_ELT(DYADVAR, 0));
	//	Rprintf("listlen =  %d\n", listlen);
	int pos = 0;
	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		double val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pConstantDyadicCovariate->value(i-1, j-1, val);
	}
	listlen = Rf_ncols(VECTOR_ELT(DYADVAR, 1));
	pos = 0;
	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		double val;
		i = missingstart[pos++];
		j = missingstart[pos++];
		val = missingstart[pos++];
		pConstantDyadicCovariate->missing(i-1, j-1, val);
	}
}

/**
 * Create one group of constant dyadic covariates
 *
 */
void setupDyadicCovariateGroup(SEXP DYADVARGROUP, Data *pData)
{
	int nDyadicCovariate = Rf_length(DYADVARGROUP);
	//    Rprintf("nDyadicCovariate %d\n", nDyadicCovariate);
	for (int dyadicCovariate = 0; dyadicCovariate < nDyadicCovariate;
			dyadicCovariate++)
	{
		SEXP as;
		PROTECT(as = Rf_install("nodeSet"));
		SEXP actorSet = Rf_getAttrib(VECTOR_ELT(DYADVARGROUP, dyadicCovariate),
				as);
		SEXP nm;
		PROTECT(nm = Rf_install("name"));
		SEXP name = Rf_getAttrib(VECTOR_ELT(DYADVARGROUP, dyadicCovariate),
				nm);
		const ActorSet * myActorSet1 = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 0)));
		const ActorSet * myActorSet2 = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 1)));
		ConstantDyadicCovariate * pConstantDyadicCovariate =
			pData->createConstantDyadicCovariate(CHAR(STRING_ELT(name, 0)),
					myActorSet1, myActorSet2);
		setupDyadicCovariate(VECTOR_ELT(DYADVARGROUP, dyadicCovariate),
				pConstantDyadicCovariate);
		SEXP mean;
		PROTECT(mean = Rf_install("mean"));
		SEXP Mean = Rf_getAttrib(VECTOR_ELT(DYADVARGROUP, dyadicCovariate),
				mean);
		pConstantDyadicCovariate->mean(REAL(Mean)[0]);
		UNPROTECT(3);
	}
}

/**
 * Unpack one set of values for a changing dyadic covariate
 *
 */
void unpackChangingDyadicPeriod(SEXP VARDYADVALS, ChangingDyadicCovariate *
		pChangingDyadicCovariate, int period)
{
	double *start = REAL(VECTOR_ELT(VARDYADVALS, 0));
	int listlen = Rf_ncols(VECTOR_ELT(VARDYADVALS, 0));
	//	Rprintf("listlen =  %d\n", listlen);
	int pos = 0;
	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		double val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pChangingDyadicCovariate->value(i - 1, j - 1, period, val);
	}
	double *missingstart = REAL(VECTOR_ELT(VARDYADVALS, 1));
	listlen = Rf_ncols(VECTOR_ELT(VARDYADVALS, 1));
	//	Rprintf("listlen =  %d\n", listlen);
	pos = 0;
	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		double val;
		i = missingstart[pos++];
		j = missingstart[pos++];
		val = missingstart[pos++];
		pChangingDyadicCovariate->missing(i - 1, j - 1, period, val);
	}
}
/**
 * Create all observations for a changing dyadic covariate
 *
 */
void setupChangingDyadicObservations(SEXP VARDYAD,
		ChangingDyadicCovariate * pChangingDyadicCovariate)
{
	int observations = Rf_length(VARDYAD);
	//   if (observations != pworkLongitudinalData->observationCount())
	// {
	//	Rf_error ("wrong number of observations in OneMode");
	//  }
	for (int period = 0; period < (observations - 1); period++)
	{
		unpackChangingDyadicPeriod(VECTOR_ELT(VARDYAD, period),
				pChangingDyadicCovariate, period);
	}
}

/**
 * Create one group of changing dyadic covariates
 *
 */
void setupChangingDyadicCovariateGroup(SEXP VARDYADGROUP, Data * pData)
{
	int nChangingDyadic = Rf_length(VARDYADGROUP);

	for (int changingDyadic = 0; changingDyadic < nChangingDyadic;
			changingDyadic++)
	{
		SEXP as;
		PROTECT(as = Rf_install("nodeSet"));
		SEXP actorSet = Rf_getAttrib(VECTOR_ELT(VARDYADGROUP, changingDyadic), as);
		SEXP nm;
		PROTECT(nm = Rf_install("name"));
		SEXP name = Rf_getAttrib(VECTOR_ELT(VARDYADGROUP, changingDyadic), nm);
		const ActorSet * myActorSet1 = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 0)));
		const ActorSet * myActorSet2 = pData->pActorSet(CHAR(STRING_ELT(
						actorSet, 1)));
		ChangingDyadicCovariate *  pChangingDyadicCovariate =
			pData->createChangingDyadicCovariate(CHAR(STRING_ELT(name, 0)),
					myActorSet1, myActorSet2);
		setupChangingDyadicObservations(VECTOR_ELT(VARDYADGROUP,
					changingDyadic),
				pChangingDyadicCovariate);
		SEXP mean;
		PROTECT(mean = Rf_install("mean"));
		SEXP Mean = Rf_getAttrib(VECTOR_ELT(VARDYADGROUP, changingDyadic),
				mean);
		pChangingDyadicCovariate->mean(REAL(Mean)[0]);
		UNPROTECT(3);
	}
}

/**
 * Create the exogenous composition change events for one actor set within
 * one group.
 */
void setupExogenousEventSet(SEXP EXOGEVENTSET, Data *pData)
{
	/* pass in the data for one actor set as two items: first
		 a list of columns: event type, period, actor, time.
		 Secondly a matrix of booleans, indicating whether active at start
		 of period.*/

	/* first find the actor set */
	SEXP as;
	PROTECT(as = Rf_install("nodeSet"));
	SEXP actorSet = Rf_getAttrib(EXOGEVENTSET, as);

	/* now process the events */
	SEXP EVENTS = VECTOR_ELT(EXOGEVENTSET, 0);
	int nEvents = Rf_length(VECTOR_ELT(EVENTS, 0));
	//Rprintf("number of rows of data frame %d\n",nEvents);
	//Rprintf("%d\n", Rf_length(EVENTS));
	int * type = INTEGER(VECTOR_ELT(EVENTS, 0));
	//Rprintf("type %d\n",*type);
	int * period = INTEGER(VECTOR_ELT(EVENTS, 1));
	//Rprintf("period %d\n",*period);
	int * actor = INTEGER(VECTOR_ELT(EVENTS, 2));
	//Rprintf("actor %d\n",*actor);
	double * time = REAL(VECTOR_ELT(EVENTS, 3));
	//Rprintf("time %5.4f\n",*time);
	const ActorSet * pActorSet = pData->pActorSet(CHAR(STRING_ELT(actorSet,
					0)));
	for (int event = 0; event < nEvents; event++)
	{
		if (*type == 1)
		{
			pData->addJoiningEvent(*period-1, pActorSet, *actor-1, *time);
		}
		else
		{
			pData->addLeavingEvent(*period-1, pActorSet, *actor-1, *time);
		}
		type++;
		period++;
		actor++;
		time++;
	}
	/* retrieve some to check*/
	//     const EventSet * myeventset= pData->pEventSet(0);
	// 	EventSet::iterator myit = myeventset->begin();
	// 	Rprintf("period 1 first event? %d %3.2f\n",(*myit)->actor(),
	// 			(*myit)->time());

	/* initialise the active flags */
	SEXP ACTIVES= VECTOR_ELT(EXOGEVENTSET, 1);

	/* this is a matrix with column for each observation and row per actor*/

	int nActors = pActorSet->n();
	int *active = LOGICAL(ACTIVES);
	for (int period = 0; period < pData->observationCount(); period++)
	{
		for (int actor = 0; actor < nActors; actor++)
		{
			pData->active(pActorSet, actor, period, *active);
			active++;
		}

	}
	UNPROTECT(1);
}

/**
 * Create one group of exogenous composition change events
 *
 */
void setupExogenousEventGroup(SEXP EXOGEVENTGROUP, Data *pData)
{
	/* pass in the data for each actor set as two items: first
		 a list of columns: event type, period, actor, time.
		 Secondly a matrix of booleans, indicating whether active at start
		 of period.*/

	int nActorSets = Rf_length(EXOGEVENTGROUP);
	//	Rprintf("number of actor sets %d\n", nActorSets);

	for (int actorSet = 0; actorSet < nActorSets; actorSet++)
	{
		setupExogenousEventSet(VECTOR_ELT(EXOGEVENTGROUP, actorSet), pData);
	}
}

/**
 *  Creates all the basic effects for one network
 */
SEXP createEffects(SEXP EFFECTS, Model *pModel, vector<Data *> * pGroupData,
		const char *networkName, int effectCol, int parmCol, int int1Col,
		int int2Col, int initValCol, int typeCol, int groupCol, int periodCol,
		int rateTypeCol, int netTypeCol, int settingCol)
{
	// find out how many effects there are
	int nEffects = Rf_length(VECTOR_ELT(EFFECTS, 0));

	// create the effects

	/* set up a vector to return the pointers in */
	SEXP effectPtrs;
	PROTECT(effectPtrs = Rf_allocVector(VECSXP, nEffects));

	for (int i = 0; i < nEffects; i++)
	{
		EffectInfo * pEffectInfo = 0;

		const char * effectName =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol), i));
		int parm1 = INTEGER(VECTOR_ELT(EFFECTS, parmCol))[i];
		double parm = parm1;
		const char * interaction1 =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int1Col), i));
		const char * interaction2 =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int2Col), i));
		double initialValue = REAL(VECTOR_ELT(EFFECTS, initValCol))[i];
		const char * effectType =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
		const char * rateType =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, rateTypeCol), i));
		const char * netType =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, netTypeCol), i));
		const char * setting =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, settingCol), i));

		if (strcmp(effectType, "rate") == 0 &&
				strcmp(effectName, "Rate") == 0)
		{
			/* find the network */
			int group = INTEGER(VECTOR_ELT(EFFECTS, groupCol))[i] - 1;
			int period = INTEGER(VECTOR_ELT(EFFECTS, periodCol))[i] - 1;

			Data * pData = (*pGroupData)[group];

			if (strcmp(setting, "") == 0)
			{
				if (!(strcmp(netType, "behavior") == 0))
				{
					NetworkLongitudinalData * pNetwork =
						pData->pNetworkData(networkName);
					pModel->basicRateParameter(pNetwork, period,
							initialValue);
				}
				else
				{
					BehaviorLongitudinalData * pNetwork =
						pData->pBehaviorData(networkName);
					pModel->basicRateParameter(pNetwork, period,
							initialValue);
				}
			}
			else
			{
				if (!(strcmp(netType, "behavior") == 0))
				{
					NetworkLongitudinalData * pNetwork =
						pData->pNetworkData(networkName);
					pModel->settingRateParameter(pNetwork, setting, period,
							initialValue);
				}
				else
				{
					Rf_error("setting found for behavior variable %s",
							networkName);
				}
			}
		}
		else if (strcmp(effectType, "rate") == 0 &&
				 strcmp(effectName, "scale") == 0)
		{
			int period = INTEGER(VECTOR_ELT(EFFECTS, periodCol))[i] - 1;
			if (strcmp(setting, "") == 0)
			{
				pModel->basicScaleParameter(period, initialValue);
			}
		else
		{
				Rf_error("setting found for variable %s", networkName);
			}
		}
		else // no rate or scale effect
		{
			pEffectInfo = pModel->addEffect(networkName,
					effectName,
					effectType,
					initialValue,
					parm,
					interaction1,
					interaction2,
					rateType);
		}

		SET_VECTOR_ELT(effectPtrs, i,
				R_MakeExternalPtr((void *) pEffectInfo,
					R_NilValue, R_NilValue));
	}

	UNPROTECT(1);
	return effectPtrs;
}

/**
 *  Creates all the interaction effects for one network
 */
SEXP createInteractionEffects(SEXP EFFECTS, Model *pModel,
		const char *networkName, int effectCol, int initValCol,
		int typeCol, int intptr1Col, int intptr2Col, int intptr3Col)
{
	// find out how many effects there are
	int nEffects = Rf_length(VECTOR_ELT(EFFECTS, 0));

	// create the effects

	/* set up a vector to return the pointers in */
	SEXP effectPtrs;
	PROTECT(effectPtrs = Rf_allocVector(VECSXP, nEffects));

	for (int i = 0; i < nEffects; i++)
	{
		EffectInfo * pEffectInfo = 0;

		const char * effectName =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol), i));
		double initialValue = REAL(VECTOR_ELT(EFFECTS, initValCol))[i];
		const char * effectType =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
		EffectInfo * pEffect1 = (EffectInfo *) R_ExternalPtrAddr(
				VECTOR_ELT(VECTOR_ELT(EFFECTS, intptr1Col), i));
		EffectInfo * pEffect2 = (EffectInfo *) R_ExternalPtrAddr(
				VECTOR_ELT(VECTOR_ELT(EFFECTS, intptr2Col), i));
		EffectInfo * pEffect3 = 0;
		if (!Rf_isNull(VECTOR_ELT(VECTOR_ELT(EFFECTS, intptr3Col), i)))
		{
			pEffect3 = (EffectInfo *) R_ExternalPtrAddr(
					VECTOR_ELT(VECTOR_ELT(EFFECTS, intptr3Col), i));
		}

		pEffectInfo = pModel->addInteractionEffect(networkName,
				effectName,
				effectType,
				initialValue,
				pEffect1,
				pEffect2,
				pEffect3);

		SET_VECTOR_ELT(effectPtrs, i,
				R_MakeExternalPtr((void *) pEffectInfo,
					R_NilValue, R_NilValue));
	}

	UNPROTECT(1);
	return effectPtrs;
}

/**
 *  Retrieves the contributions to all possible tie flips or behavior changes for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getChangeContributionStatistics(SEXP EFFECTSLIST,
		const StatisticCalculator * pCalculator, vector<vector<double *> > *rChangeContributions)
{

	// get the column names from the names attribute
	SEXP cols = PROTECT(Rf_install("names"));
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

	for (int ii = 0; ii < Rf_length(EFFECTSLIST); ii++)
	{
		SEXP EFFECTS = VECTOR_ELT(EFFECTSLIST, ii);

		for (int i = 0; i < Rf_length(VECTOR_ELT(EFFECTS,0)); i++)
		{
			const char * effectType = CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
			const char * netType = CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, netTypeCol), i));
			if(strcmp(netType, "oneMode") == 0 || strcmp(netType, "bipartite") == 0 ||
									strcmp(netType, "behavior") == 0)
			{
				// todo At the moment, change contributions cannot be calculated for endowment or creation effects
				// modifications in the corresponding methods (calculateNetworkEndowmentStatistics, calculateNetworkCreationStatistics,
				// and calculateBehaviorStatistics) in StatisticCalculator.cpp would be necessary!!!
				//if (strcmp(effectType, "eval") == 0 || strcmp(effectType, "endow") == 0 || strcmp(effectType, "creation") == 0)
				if (strcmp(effectType, "eval") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *) R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(EFFECTS,pointerCol), i));
					if(rChangeContributions != 0)
					{
						rChangeContributions->push_back(pCalculator->staticChangeContributions(pEffectInfo));
					}
				}
			}
		}
	}
	UNPROTECT(1);
}

/**
 *  Retrieves the statistics of individual actors for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getActorStatistics(SEXP EFFECTSLIST,
		const StatisticCalculator * pCalculator, vector<double *> *rActorStatistics)
{

	// get the column names from the names attribute
	SEXP cols = PROTECT(Rf_install("names"));
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

	for (int ii = 0; ii < Rf_length(EFFECTSLIST); ii++)
	{
		SEXP EFFECTS = VECTOR_ELT(EFFECTSLIST, ii);

		for (int i = 0; i < Rf_length(VECTOR_ELT(EFFECTS,0)); i++)
		{
			const char * effectType = CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
			const char * netType = CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, netTypeCol), i));
			if(strcmp(netType, "oneMode") == 0 || strcmp(netType, "behavior") == 0)
			{
				if (strcmp(effectType, "eval") == 0 || strcmp(effectType, "endow") == 0 || strcmp(effectType, "creation") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *) R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(EFFECTS,pointerCol), i));
					if(rActorStatistics != 0)
					{
						rActorStatistics->push_back(pCalculator->actorStatistics(pEffectInfo));
					}
				}
			}
		}
	}
	UNPROTECT(1);
}

/**
 *  Retrieves the values of the statistics and scores for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getStatistics(SEXP EFFECTSLIST,
		const StatisticCalculator * pCalculator,
		int period, int group, const Data *pData,
		const EpochSimulation * pEpochSimulation,
		vector<double> * rfra, vector<double> *rscore)
{

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


	double statistic = 0;
	double score = 0;
	int istore = 0;

	for (int ii = 0; ii < Rf_length(EFFECTSLIST); ii++)
	{
		const char * networkName =
			CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, ii),
							nameCol), 0));
		SEXP EFFECTS = VECTOR_ELT(EFFECTSLIST, ii);

		for (int i = 0; i < Rf_length(VECTOR_ELT(EFFECTS,0)); i++)
		{
			const char * effectName =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol),  i));
			//	int parm = INTEGER(VECTOR_ELT(EFFECTS, parmCol))[i];
			const char * interaction1 =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int1Col),i));
			//	const char * interaction2 =
			//CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int2Col), i));
			//double initialValue =
			//REAL(VECTOR_ELT(EFFECTS, initValCol))[i];
			const char * effectType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
			const char * netType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, netTypeCol), i));
			const char * rateType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, rateTypeCol), i));
			const char * setting =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, settingCol), i));
			//	Rprintf("%s %s \n", effectType, netType);
			if (strcmp(effectType, "rate") == 0)
			{
				if (strcmp(effectName, "Rate") == 0 ||
					strcmp(effectName, "scale") == 0)
				{
					int groupno =
						INTEGER(VECTOR_ELT(EFFECTS, groupCol))[i];
					int periodno =
						INTEGER(VECTOR_ELT(EFFECTS, periodCol))[i];
					if ((periodno - 1) == period && (groupno - 1) == group)
					{

						if (strcmp(netType, "behavior") == 0)
						{
							LongitudinalData * pNetworkData =
								pData->pBehaviorData(networkName);
							statistic = pCalculator->distance(pNetworkData,
									period);
							//	Rprintf("%f behavior dist\n", statistic);
							if (pEpochSimulation)
							{
								const DependentVariable * pVariable =
									pEpochSimulation->pVariable(networkName);
								score = pVariable->basicRateScore();
							}
							else
							{
								score = 0;
							}
						}
						else if (strcmp(netType, "continuous") == 0)
						{
							statistic = pCalculator->totalDistance(period);

							if (pEpochSimulation)
							{
								score = pEpochSimulation->pSdeSimulation()->basicScaleScore();
								// Rprintf("The tau score gets update with %f\n", score); NYNKE
							}
							else
							{
								score = 0;
							}
						}
						else
						{
							if (strcmp(setting, "") == 0)
							{
								LongitudinalData * pNetworkData =
									pData->pNetworkData(networkName);
								statistic = pCalculator->distance(pNetworkData,
										period);
								if (pEpochSimulation)
								{
									/* find  dependent variable for the score */
									const DependentVariable * pVariable =
										pEpochSimulation->
										pVariable(networkName);
									score = pVariable->basicRateScore();
								}
								else
								{
									score = 0;
								}
							}
							else
							{
								LongitudinalData * pNetworkData =
									pData->pNetworkData(networkName);
								statistic =
									pCalculator->settingDistance(pNetworkData,
											setting, period);
								if (pEpochSimulation)
								{
									/* find dependent variable for the score */
									const DependentVariable * pVariable =
										pEpochSimulation->
										pVariable(networkName);
									score =
										pVariable->settingRateScore(setting);
								}
								else
								{
									score = 0;
								}

							}
						}
					}
					else
					{
						statistic = 0;
						score = 0;
					}
				}
				else if (strcmp(rateType, "structural") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
								VECTOR_ELT(VECTOR_ELT(EFFECTS,
										pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);
					if (pEpochSimulation)
					{
						const DependentVariable * pVariable =
							pEpochSimulation->pVariable(networkName);
						const NetworkVariable * pNetworkVariable;
						if (strcmp(interaction1, "") == 0)
						{
							pNetworkVariable =
								(const NetworkVariable *)
								pEpochSimulation->pVariable(networkName);
						}
						else
						{
							pNetworkVariable =
								(const NetworkVariable *)
								pEpochSimulation->pVariable(interaction1);
						}
						if (strcmp(effectName, "outRate") == 0)
						{
							score =
								pVariable->outDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "inRate") == 0)
						{
							score =
								pVariable->inDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "recipRate") == 0)
						{
							score =
								pVariable->reciprocalDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "outRateInv") == 0)
						{
							score =
								pVariable->inverseOutDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "outRateLog") == 0)
						{
							score =
								pVariable->logOutDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "inRateInv") == 0)
						{
							score =
								pVariable->inverseInDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "inRateLog") == 0)
						{
							score =
								pVariable->logInDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "outRateLog") == 0)
						{
							score =
								pVariable->logOutDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "recipRateInv") == 0)
						{
							score =
								pVariable->inversereciprocalDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "recipRateLog") == 0)
						{
							score =
								pVariable->logreciprocalDegreeScore(pNetworkVariable);
						}
						else
						{

							Rf_error("Unexpected rate effect %s\n",
									effectName);
						}
					}
					else
					{
						score = 0;
					}
				}
				else if (strcmp(rateType, "diffusion") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
								VECTOR_ELT(VECTOR_ELT(EFFECTS,
										pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);
					if (pEpochSimulation)
					{
						if (strcmp(effectName, "avExposure") == 0 ||
								strcmp(effectName, "totExposure") == 0 ||
								strcmp(effectName, "susceptAvIn") == 0 ||
								strcmp(effectName, "infectIn") == 0 ||
								strcmp(effectName, "infectDeg") == 0 ||
								strcmp(effectName, "infectOut") == 0 ||
								strcmp(effectName, "susceptAvCovar") == 0 ||
								strcmp(effectName, "infectCovar") == 0)
						{
							score = pEpochSimulation->score(pEffectInfo);
						}
						else
						{
							Rf_error("Unexpected rate effect %s\n",
									effectName);
						}
					}
					else
					{
						score = 0;
					}
				}
				else
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
								VECTOR_ELT(VECTOR_ELT(EFFECTS,
										pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);

					if (pEpochSimulation)
					{
						ConstantCovariate * pConstantCovariate =
							pData->pConstantCovariate(interaction1);
						ChangingCovariate * pChangingCovariate =
							pData->pChangingCovariate(interaction1);
						BehaviorVariable * pBehavior =
							(BehaviorVariable *)
							pEpochSimulation->pVariable(interaction1);
						//find the network

						const DependentVariable * pVariable =
							pEpochSimulation->pVariable(networkName);

						if (pConstantCovariate)
						{
							score = pVariable->constantCovariateScore(
									pConstantCovariate);
						}
						else if (pChangingCovariate)
						{
							score = pVariable->changingCovariateScore(
									pChangingCovariate);
						}
						else if (pBehavior)
						{
							score = pVariable->behaviorVariableScore(
									pBehavior);
						}
						else
						{
							Rf_error("No individual covariate named %s.",
									interaction1);
						}
					}
					else
					{
						score = 0;
					}
				}

			}
			else
			{
				if (strcmp(effectType, "eval") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
								VECTOR_ELT(VECTOR_ELT(EFFECTS,
										pointerCol), i));
					statistic =	pCalculator->statistic(pEffectInfo);
					if (pEpochSimulation)
					{
						score = pEpochSimulation->score(pEffectInfo);
					}
					else
					{
						score = 0;
					}
				}
				else if (strcmp(effectType, "endow") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
								VECTOR_ELT(VECTOR_ELT(EFFECTS,
										pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);
					if (strcmp(netType, "behavior") != 0)
					{
						statistic = -1 * statistic;
					}
					if (pEpochSimulation)
					{
						score = pEpochSimulation->score(pEffectInfo);
					}
					else
					{
						score = 0;
					}
				}
				else if (strcmp(effectType, "creation") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
								VECTOR_ELT(VECTOR_ELT(EFFECTS,
										pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);
					if (pEpochSimulation)
					{
						score = pEpochSimulation->score(pEffectInfo);
					}
					else
					{
						score = 0;
					}
				}
				else if (strcmp(effectType, "gmm") == 0)
				{
				  EffectInfo * pEffectInfo = (EffectInfo *)
				  R_ExternalPtrAddr(
				    VECTOR_ELT(VECTOR_ELT(EFFECTS,
                              pointerCol), i));
				  statistic = pCalculator->statistic(pEffectInfo);
				  if (pEpochSimulation)
				  {
				    score = pEpochSimulation->score(pEffectInfo);
				  }
				  else
				  {
				    score = 0;
				  }
				}
				else
				{
					Rf_error("invalid effect type %s\n", effectType);
				}
			}
			(*rfra)[istore] = statistic;
			if (pEpochSimulation)
			{
				//		Rprintf("%f %d \n", score, istore);
				(*rscore)[istore] = score;
			}
			istore++; /* keep forgetting to move the ++ */
		}
	}
	UNPROTECT(1);
}

/**
 *  retrieves the values of the scores and derivatives for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Only used in maximum likelihood.
 */
void getScores(SEXP EFFECTSLIST, int period, int group,
		const MLSimulation * pMLSimulation,
		vector<double> * rderiv, vector<double> *rscore)
{

	// get the column names from the names attribute
	SEXP cols = PROTECT(Rf_install("names"));
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

	int storescore = 0;
	int storederiv = 0;

	for (int ii = 0; ii < Rf_length(EFFECTSLIST); ii++)
	{
		const char * networkName =
			CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, ii),
							nameCol), 0));
		SEXP EFFECTS = VECTOR_ELT(EFFECTSLIST, ii);

		for (int i = 0; i < Rf_length(VECTOR_ELT(EFFECTS,0)); i++)
		{
			const char * effectName =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol),  i));
			const char * effectType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
			if (strcmp(effectType, "rate") == 0)
			{
				if (strcmp(effectName, "Rate") == 0)
				{
					int groupno =
						INTEGER(VECTOR_ELT(EFFECTS, groupCol))[i];
					int periodno =
						INTEGER(VECTOR_ELT(EFFECTS, periodCol))[i];
					if ((periodno - 1) == period && (groupno - 1) == group)
					{
						const DependentVariable * pVariable =
							pMLSimulation->pVariable(networkName);
						(*rscore)[storescore++] = pVariable->basicRateScore();
						(*rderiv)[storederiv++] =
							pVariable->basicRateDerivative();
					}
					else
					{
						(*rscore)[storescore++] = 0;
						(*rderiv)[storederiv++] = 0;
					}
				}
				else
				{
					Rf_error("Non constant rate effects are not yet %s",
							"implemented for maximum likelihood.");
				}
			}
			else
			{
				EffectInfo * pEffectInfo = (EffectInfo *)
					R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(EFFECTS,
									pointerCol), i));
				(*rscore)[storescore++] = pMLSimulation->score(pEffectInfo);

				// get the map of derivatives
				map<const EffectInfo *, double > deriv =
					pMLSimulation->derivative(pEffectInfo);

				for (int j = 0; j < Rf_length(VECTOR_ELT(EFFECTS,0)); j++)
				{
					const char * effectType =
						CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), j));

					if (!(strcmp(effectType, "rate") == 0))
					{
						//	Rprintf("%s %s \n", effectType, netType);
						EffectInfo * pEffectInfo2 = (EffectInfo *)
							R_ExternalPtrAddr(
									VECTOR_ELT(VECTOR_ELT(EFFECTS,
											pointerCol), j));

						(*rderiv)[storederiv++] =
							pMLSimulation->derivative(pEffectInfo,
									pEffectInfo2);
					}
				}
			}
		}
	}
	UNPROTECT(1);
}


