/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07utilities.cpp
 *
 * Description: This module contains various utilities, including:
 *
 *  1) printOutData: dump out the data from a single data object for
 *  profiling in sienaProfile.exe. Called from setupModelOptions if
 *  requested by user.
 *  NB. This function is probably not up to date or complete.
 *
 *  2) getMiniStepDF: create a data frame as a SEXP from a ministep
 *
 *  3) getChainDF: create a data frame as a SEXP from a chain
 *
 *  4) getMiniStepList: create a list format SEXP from a ministep
 *
 *  5) getChainList: create a list format SEXP from a chain
 *
 *  SEXP's can be printed within C using PrintValue(SEXP x)
 *****************************************************************************/
/**
 * @file
 * Dump data out for profiling.
 * Convert ministeps and chains to R format objects.
 */
#include <stdexcept>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>

#include "siena07utilities.h"

#include "siena07internals.h"
#include "data/Data.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "data/LongitudinalData.h"
#include "data/NetworkLongitudinalData.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/ConstantCovariate.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "utils/Utils.h"
#include "model/EpochSimulation.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/Chain.h"
#include "model/ml/MLSimulation.h"
#include "model/ml/MiniStep.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/BehaviorChange.h"
#include "model/State.h"

using namespace std;
using namespace siena;

//--------------------------------------
// utility functions to process groups
//-------------------------------------

/** Calculate the period number of this group and period, to access
 * correct stored chain.
 */
int periodFromStart(vector<Data *> & pGroupData, int group, int period)
{

	int periodFromStart = 0;
	for (int i = 0; i < group; i++)
	{
		periodFromStart += pGroupData[i]->observationCount() - 1;
	}
	periodFromStart += period;

	return periodFromStart;
}

/** Calculate the total number of periods in all groups, which is the dimension
 * of some returned arrays.
 */

int totalPeriods(vector<Data *> & pGroupData)
{
	int nGroups = pGroupData.size();

	int totObservations = 0;

	for (int group = 0; group < nGroups; group++)
	{
		totObservations += pGroupData[group]->observationCount() - 1;
	}
	return totObservations;

}

/**
 * Traps errors so R can stop the function rather than being stoppped itself.
 *
 */
void Rterminate()
{
	try
	{
		throw;
	}
	catch(exception& e)
	{
		error(e.what());
	}
}


/**
 * print out the data for profiling with gprof
 *
 */
void printOutData(Data *pData)
{

	ofstream myfile ("data.txt");

	if (myfile.is_open())
	{
		myfile << pData->observationCount() << endl;
		const std::vector<LongitudinalData * > rVariables =
			pData->rDependentVariableData();
		int nActors = rVariables[0]->n();

		myfile << rVariables[0]->n();
		for (unsigned i = 0; i < pData->rDependentVariableData().size(); i++)
		{
			NetworkLongitudinalData * pNetworkData =
				dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
			BehaviorLongitudinalData * pBehaviorData =
				dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);

			myfile << rVariables[i]->name() << endl;

			for(int period = 0; period < pData->observationCount(); period++)
			{
				if (pNetworkData)
				{
					if (period == 0) myfile << "oneMode" << endl;
					myfile << pNetworkData->pNetwork(period)->tieCount() << endl;
					for (TieIterator
							 iter=pNetworkData-> pNetwork(period)->ties();
						 iter.valid();
						 iter.next())
					{
						myfile << iter.ego() << " "
							   << iter.alter() << " "
						 << iter.value() << endl;
					}
					myfile << pNetworkData->
						pMissingTieNetwork(period)->tieCount() << endl;
					for (TieIterator
							 iter=pNetworkData->
							 pMissingTieNetwork(period)->ties();
						 iter.valid();
						 iter.next())
					{
						myfile << iter.ego() << " "
							   << iter.alter() << " "
						 << iter.value() << endl;
					}
					myfile << pNetworkData->
						pStructuralTieNetwork(period)->tieCount() << endl;
					for (TieIterator
							iter=pNetworkData->
							pStructuralTieNetwork(period)->ties();
							iter.valid();
							iter.next())
					{
						myfile << iter.ego() << " "
							<< iter.alter() << " "
							<< iter.value() << endl;
					}
					// other attributes uponly downonly
					myfile << pNetworkData->upOnly(period) <<  " " <<
						pNetworkData->downOnly(period) << endl;
				}

				else if (pBehaviorData)
				{
					if (period ==0 )
					{
						myfile << "behavior" << endl;
						myfile <<  pBehaviorData->n() << endl;
					}
					for (int ii = 0; ii < pBehaviorData->n(); ii++)
					{
						myfile << ii << " " <<
							pBehaviorData->value(period, ii) << endl;
					}
					for (int ii = 0; ii < pBehaviorData->n(); ii++)
					{
						myfile << ii << " " <<
							pBehaviorData->missing(period, ii) << endl;
					}
					// other attributes similarityMean uponly downonly
					myfile << pBehaviorData->similarityMean() << endl;
					myfile << pBehaviorData->upOnly(period) << " " <<
						pBehaviorData->downOnly(period) << endl;
				}
				else
				{
					throw ("Unexpected class of dependent variable");
				}

			}
		}
		const std::vector<ConstantCovariate * > rConstantCovariates =
			pData->rConstantCovariates();
		for (unsigned i = 0; i < pData->rConstantCovariates().size(); i++)
		{
			ConstantCovariate * pCovariate = rConstantCovariates[i];

			myfile << pCovariate->name() << endl;
			myfile << "constantcovariate" << endl;
			myfile <<  nActors << endl;

			for (int ii = 0; ii < nActors; ii++)
			{
				myfile << ii << " " <<	pCovariate->value(ii) << endl;
			}
			for (int ii = 0; ii < nActors; ii++)
			{
				myfile << ii << " " <<	pCovariate->missing(ii) << endl;
			}

			myfile << pCovariate->similarityMean() << endl;
			myfile << pCovariate->range() << endl;
		}
		const std::vector<ChangingCovariate * > rChangingCovariates =
			pData->rChangingCovariates();

		for (unsigned i = 0; i < pData->rChangingCovariates().size(); i++)
		{
			ChangingCovariate * pCovariate = rChangingCovariates[i];

			myfile << pCovariate->name() << endl;
			myfile << "changingcovariate" << endl;
			myfile <<  nActors << endl;

			for (int period = 0; period < pData->observationCount() - 1;
				 period++)
			{
				for (int ii = 0; ii < nActors; ii++)
				{
					myfile << ii << " " <<
						pCovariate->value(ii, period) << endl;
				}
				for (int ii = 0; ii < nActors; ii++)
				{
					myfile << ii << " " <<
						pCovariate->missing(ii, period) << endl;
				}
			}
			myfile << pCovariate->similarityMean() << endl;
			myfile << pCovariate->range() << endl;
		}

		const std::vector<ConstantDyadicCovariate * >
			rConstantDyadicCovariates =	pData->rConstantDyadicCovariates();

		for (unsigned i = 0; i < pData->rConstantDyadicCovariates().size(); i++)
		{
			ConstantDyadicCovariate * pCovariate = rConstantDyadicCovariates[i];

			myfile << pCovariate->name() << endl;
			myfile << "constantdyadiccovariate" << endl;
			myfile <<  nActors << endl;

			for (int j = 0; j < nActors; j++)
			{
				for (int k = 0; k < nActors; k++)
				{
					myfile << j << " " << k << " " << pCovariate->value(j, k)
						   << endl;
					myfile << j << " " << k << " " <<
						pCovariate->missing(j, k) << endl;
				}
			}
			myfile << pCovariate->mean() << endl;
		}
		const std::vector<ChangingDyadicCovariate * >
			rChangingDyadicCovariates =	pData->rChangingDyadicCovariates();

		for (unsigned i = 0; i < pData->rChangingDyadicCovariates().size(); i++)
		{
			ChangingDyadicCovariate * pCovariate = rChangingDyadicCovariates[i];

			myfile << pCovariate->name() << endl;
			myfile << "changingdyadiccovariate" << endl;
			myfile <<  nActors << endl;

			for (int period = 0; period < pData->observationCount() - 1;
				 period++)
			{
				for (int j = 0; j < nActors; j++)
				{
					for (int k = 0; k < nActors; k++)
					{
						myfile << j << " " << k << " " <<
							pCovariate->value(j, k, period)
							   << endl;
						myfile << j << " " << k << " " <<
							pCovariate->missing(j, k, period) << endl;
					}
				}
			}
			myfile << pCovariate->mean() << endl;
		}
	}
}

/** Create an R vector from a behavior variable for a single period
 *
 */
SEXP getBehaviorValues(const BehaviorVariable & behavior)
{
    int n = behavior.n();
	SEXP ans = PROTECT(allocVector(INTSXP, n));
    int *ians = INTEGER(ans);
	const int *pValues = behavior.values();
    for (int i = 0; i < n; i++)
	{
		ians[i] = pValues[i];
	}
	UNPROTECT(1);
    return ans  ;
}

/**
 * Create an R matrix from a network variable for a single period
 */
SEXP getAdjacency(const Network& net)
{
    int n=net.n();
    int m=net.m();
	SEXP ans = PROTECT(allocMatrix(INTSXP, n, m));
    int *ians = INTEGER(ans);
    /* initialise the memory: possibly only necessary in case of error! */
    for (int i = 0; i<n*m;i++)
	ians[i]=0;
    for (TieIterator iter=net.ties(); iter.valid(); iter.next())
    {
		ians[iter.ego()+(iter.alter()*n)] = iter.value();
    }

    UNPROTECT(1);
    return ans;
}

/**
 * Create an R 3 column matrix from a network variable for a single period
 */
SEXP getEdgeList(const Network& net)
{
	int nties = net.tieCount();
	SEXP ans = PROTECT(allocMatrix(INTSXP, nties, 3));
    int *ians = INTEGER(ans);
    /* initialise the memory: possibly only neccesary in case of error! */
	for (int i = 0; i < nties * 3; i++) {
		ians[i] = 0;
	}
	int irow = 0;
	for (TieIterator iter = net.ties(); iter.valid(); iter.next()) {
		ians[irow ] = iter.ego() + 1;
		ians[nties + irow] = iter.alter() + 1;
		ians[2 * nties + irow] = iter.value();
		irow ++;
    }

    UNPROTECT(1);
    return ans;
}
/**
 * utilities to access chains and ministeps
 */
namespace siena
{

SEXP net_to_sexp(const Network * pNet) {
	return getEdgeList(*pNet);
}

SEXP var_to_sexp(DependentVariable * pVar) {
	NetworkVariable * pNetworkVariable = dynamic_cast<NetworkVariable *>(pVar);
	BehaviorVariable * pBehaviorVariable = dynamic_cast<BehaviorVariable *>(pVar);
	if (pNetworkVariable) {
		return net_to_sexp(pNetworkVariable->pNetwork());
	} else if (pBehaviorVariable) {
		return getBehaviorValues(*pBehaviorVariable);
	} else {
		throw domain_error(pVar->name() + ": unexpected class of DependentVariable");
	}
}

/** Create a data frame with a single row from a ministep. (prints nicely
 * with PrintValue)
 */
SEXP getMiniStepDF(const MiniStep& miniStep)
{
	SEXP MINISTEP, classname, dimnames, colnames;
	if (miniStep.networkMiniStep() || miniStep.behaviorMiniStep())
	{
		PROTECT(colnames = allocVector(STRSXP, 10));
		SET_STRING_ELT(colnames, 0, mkChar("Aspect"));
		SET_STRING_ELT(colnames, 1, mkChar("Var"));
		SET_STRING_ELT(colnames, 2, mkChar("VarName"));
		SET_STRING_ELT(colnames, 3, mkChar("Ego"));
		SET_STRING_ELT(colnames, 4, mkChar("Alter"));
		SET_STRING_ELT(colnames, 5, mkChar("Diff"));
		SET_STRING_ELT(colnames, 6, mkChar("ReciRate"));
		SET_STRING_ELT(colnames, 7, mkChar("LogOptionSetProb"));
		SET_STRING_ELT(colnames, 8, mkChar("LogChoiceProb"));
		SET_STRING_ELT(colnames, 9, mkChar("Diagonal"));

		PROTECT(MINISTEP = allocVector(VECSXP, 10));

		if (miniStep.networkMiniStep())
		{
			const NetworkChange& networkChange =
				dynamic_cast<const NetworkChange &>(miniStep);
			SET_VECTOR_ELT(MINISTEP, 0, mkString("Network"));
			SET_VECTOR_ELT(MINISTEP, 4, ScalarInteger(networkChange.alter()));
			SET_VECTOR_ELT(MINISTEP, 5, ScalarInteger(0));
		}
		else
		{
			const BehaviorChange& behaviorChange =
				dynamic_cast<const BehaviorChange &>(miniStep);
			SET_VECTOR_ELT(MINISTEP, 0, mkString("Behavior"));
			SET_VECTOR_ELT(MINISTEP, 4, ScalarInteger(0));
			SET_VECTOR_ELT(MINISTEP, 5,
				ScalarInteger(behaviorChange.difference()));
		}
		SET_VECTOR_ELT(MINISTEP, 1, ScalarInteger(miniStep.variableId()));
		SET_VECTOR_ELT(MINISTEP, 2, mkString(miniStep.variableName().c_str()));
		SET_VECTOR_ELT(MINISTEP, 3, ScalarInteger(miniStep.ego()));
		SET_VECTOR_ELT(MINISTEP, 6, ScalarReal(miniStep.reciprocalRate()));
		SET_VECTOR_ELT(MINISTEP, 7,
			ScalarReal(miniStep.logOptionSetProbability()));
		SET_VECTOR_ELT(MINISTEP, 8,
			ScalarReal(miniStep.logChoiceProbability()));
		SET_VECTOR_ELT(MINISTEP, 9,
			ScalarLogical(miniStep.diagonal()));

		namesgets(MINISTEP, colnames);

		PROTECT(dimnames = allocVector(INTSXP, 2));
		int * idimnames = INTEGER(dimnames);
		idimnames[0] = NA_INTEGER;
		idimnames[1] = -1;
		setAttrib(MINISTEP, R_RowNamesSymbol, dimnames);

		PROTECT(classname = allocVector(STRSXP, 1));
		SET_STRING_ELT(classname, 0, mkChar("data.frame"));
		classgets(MINISTEP, classname);

		UNPROTECT(4);
		return MINISTEP;
	}
	else
		return R_NilValue;
}
/**
 * Create a data frame from a chain. (prints nicely with PrintValue)
 */
SEXP getChainDF(const Chain& chain, bool sort)
{
	SEXP ans, col0, col1, col2, col3, col4, col5, col6, col7, col8, col9,
		colnames, dimnames, classname;
	PROTECT(colnames = allocVector(STRSXP, 10));
	SET_STRING_ELT(colnames, 0, mkChar("Aspect"));
	SET_STRING_ELT(colnames, 1, mkChar("Var"));
	SET_STRING_ELT(colnames, 2, mkChar("VarName"));
	SET_STRING_ELT(colnames, 3, mkChar("Ego"));
	SET_STRING_ELT(colnames, 4, mkChar("Alter"));
	SET_STRING_ELT(colnames, 5, mkChar("Diff"));
	SET_STRING_ELT(colnames, 6, mkChar("ReciRate"));
	SET_STRING_ELT(colnames, 7, mkChar("LogOptionSetProb"));
	SET_STRING_ELT(colnames, 8, mkChar("LogChoiceProb"));
	SET_STRING_ELT(colnames, 9, mkChar("Diagonal"));

	PROTECT(ans = allocVector(VECSXP, 10));
	int numberRows = chain.ministepCount() - 1;
	PROTECT(col0 = allocVector(STRSXP, numberRows));

	PROTECT(col1 = allocVector(INTSXP, numberRows));
	int * icol1 = INTEGER(col1);
	PROTECT(col2 = allocVector(STRSXP, numberRows));
	PROTECT(col3 = allocVector(INTSXP, numberRows));
	int * icol3 = INTEGER(col3);
	PROTECT(col4 = allocVector(INTSXP, numberRows));
	int * icol4 = INTEGER(col4);
	PROTECT(col5 = allocVector(INTSXP, numberRows));
	int * icol5 = INTEGER(col5);
	PROTECT(col6 = allocVector(REALSXP, numberRows));
	double * rcol6 = REAL(col6);
	PROTECT(col7 = allocVector(REALSXP, numberRows));
	double * rcol7 = REAL(col7);
	PROTECT(col8 = allocVector(REALSXP, numberRows));
	double * rcol8 = REAL(col8);
	PROTECT(col9 = allocVector(LGLSXP, numberRows));
	int * icol9 = LOGICAL(col9);

	MiniStep *pMiniStep = chain.pFirst()->pNext();
	for (int i = 0; i < numberRows; i++)
	{
		SEXP ministep;
		PROTECT(ministep = getMiniStepDF(*pMiniStep));
		//put them in the data frame
		//	PrintValue(VECTOR_ELT(ministep, 0));
		SET_STRING_ELT(col0, i, STRING_ELT(VECTOR_ELT(ministep, 0), 0));
		icol1[i] =  INTEGER(VECTOR_ELT(ministep, 1))[0];
		SET_STRING_ELT(col2, i, STRING_ELT(VECTOR_ELT(ministep, 2), 0));
		icol3[i] =  INTEGER(VECTOR_ELT(ministep, 3))[0];
		icol4[i] =  INTEGER(VECTOR_ELT(ministep, 4))[0];
		icol5[i] =  INTEGER(VECTOR_ELT(ministep, 5))[0];
		rcol6[i] =  REAL(VECTOR_ELT(ministep, 6))[0];
		rcol7[i] =  REAL(VECTOR_ELT(ministep, 7))[0];
		rcol8[i] =  REAL(VECTOR_ELT(ministep, 8))[0];
		icol9[i] =  LOGICAL(VECTOR_ELT(ministep, 9))[0];
		pMiniStep = pMiniStep->pNext();
		UNPROTECT(1);
	}
	SET_VECTOR_ELT(ans, 0, col0);
	SET_VECTOR_ELT(ans, 1, col1);
	SET_VECTOR_ELT(ans, 2, col2);
	SET_VECTOR_ELT(ans, 3, col3);
	SET_VECTOR_ELT(ans, 4, col4);
	SET_VECTOR_ELT(ans, 5, col5);
	SET_VECTOR_ELT(ans, 6, col6);
	SET_VECTOR_ELT(ans, 7, col7);
	SET_VECTOR_ELT(ans, 8, col8);
	SET_VECTOR_ELT(ans, 9, col9);

	namesgets(ans, colnames);

	PROTECT(dimnames = allocVector(INTSXP, 2));
	int * idimnames = INTEGER(dimnames);
	idimnames[0] = NA_INTEGER;
	idimnames[1] = -numberRows;
	setAttrib(ans, R_RowNamesSymbol, dimnames);

	PROTECT(classname = allocVector(STRSXP, 1));
	SET_STRING_ELT(classname, 0, mkChar("data.frame"));
	classgets(ans, classname);

	// try to sort it by variable, ego and alter
	SEXP R_fcall1, ordering, R_fcall2, ansnew;
	PROTECT(R_fcall1 = lang4(install("order"), col1, col3, col4));
	PROTECT(ordering = eval(R_fcall1, R_GlobalEnv));
	// now sort the data frame using [.data.frame and ordering
	PROTECT(R_fcall2 = lang4(install("[.data.frame"),
			ans, ordering, R_MissingArg));
	PROTECT(ansnew = eval(R_fcall2, R_GlobalEnv));

	UNPROTECT(18);
	if (sort)
	{
		return ansnew;
	}
	else
	{
		return ans;
	}
}
/**
 * Create a data frame from a vector of ministeps.
 */
SEXP getDFFromVector(const vector< MiniStep *>& rMiniSteps, bool sort)
{
	SEXP ans, col0, col1, col2, col3, col4, col5, col6, col7, col8, col9,
		colnames, dimnames, classname;
	PROTECT(colnames = allocVector(STRSXP, 10));
	SET_STRING_ELT(colnames, 0, mkChar("Aspect"));
	SET_STRING_ELT(colnames, 1, mkChar("Var"));
	SET_STRING_ELT(colnames, 2, mkChar("VarName"));
	SET_STRING_ELT(colnames, 3, mkChar("Ego"));
	SET_STRING_ELT(colnames, 4, mkChar("Alter"));
	SET_STRING_ELT(colnames, 5, mkChar("Diff"));
	SET_STRING_ELT(colnames, 6, mkChar("ReciRate"));
	SET_STRING_ELT(colnames, 7, mkChar("LogOptionSetProb"));
	SET_STRING_ELT(colnames, 8, mkChar("LogChoiceProb"));
	SET_STRING_ELT(colnames, 9, mkChar("Diagonal"));

	PROTECT(ans = allocVector(VECSXP, 10));
	int numberRows = rMiniSteps.size();
	PROTECT(col0 = allocVector(STRSXP, numberRows));

	PROTECT(col1 = allocVector(INTSXP, numberRows));
	int * icol1 = INTEGER(col1);
	PROTECT(col2 = allocVector(STRSXP, numberRows));
	PROTECT(col3 = allocVector(INTSXP, numberRows));
	int * icol3 = INTEGER(col3);
	PROTECT(col4 = allocVector(INTSXP, numberRows));
	int * icol4 = INTEGER(col4);
	PROTECT(col5 = allocVector(INTSXP, numberRows));
	int * icol5 = INTEGER(col5);
	PROTECT(col6 = allocVector(REALSXP, numberRows));
	double * rcol6 = REAL(col6);
	PROTECT(col7 = allocVector(REALSXP, numberRows));
	double * rcol7 = REAL(col7);
	PROTECT(col8 = allocVector(REALSXP, numberRows));
	double * rcol8 = REAL(col8);
	PROTECT(col9 = allocVector(LGLSXP, numberRows));
	int * icol9 = INTEGER(col9);

	for (int i = 0; i < numberRows; i++)
	{
		SEXP ministep;
		PROTECT(ministep = getMiniStepDF(*rMiniSteps[i]));
		//put them in the data frame
		//	PrintValue(VECTOR_ELT(ministep, 0));
		SET_STRING_ELT(col0, i, STRING_ELT(VECTOR_ELT(ministep, 0), 0));
		icol1[i] =  INTEGER(VECTOR_ELT(ministep, 1))[0];
		SET_STRING_ELT(col2, i, STRING_ELT(VECTOR_ELT(ministep, 2), 0));
		icol3[i] =  INTEGER(VECTOR_ELT(ministep, 3))[0];
		icol4[i] =  INTEGER(VECTOR_ELT(ministep, 4))[0];
		icol5[i] =  INTEGER(VECTOR_ELT(ministep, 5))[0];
		rcol6[i] =  REAL(VECTOR_ELT(ministep, 6))[0];
		rcol7[i] =  REAL(VECTOR_ELT(ministep, 7))[0];
		rcol8[i] =  REAL(VECTOR_ELT(ministep, 8))[0];
		icol9[i] =  INTEGER(VECTOR_ELT(ministep, 9))[0];
		UNPROTECT(1);
	}
	SET_VECTOR_ELT(ans, 0, col0);
	SET_VECTOR_ELT(ans, 1, col1);
	SET_VECTOR_ELT(ans, 2, col2);
	SET_VECTOR_ELT(ans, 3, col3);
	SET_VECTOR_ELT(ans, 4, col4);
	SET_VECTOR_ELT(ans, 5, col5);
	SET_VECTOR_ELT(ans, 6, col6);
	SET_VECTOR_ELT(ans, 7, col7);
	SET_VECTOR_ELT(ans, 8, col8);
	SET_VECTOR_ELT(ans, 9, col9);

	namesgets(ans, colnames);

	PROTECT(dimnames = allocVector(INTSXP, 2));
	int * idimnames = INTEGER(dimnames);
	idimnames[0] = NA_INTEGER;
	idimnames[1] = -numberRows;
	setAttrib(ans, R_RowNamesSymbol, dimnames);

	PROTECT(classname = allocVector(STRSXP, 1));
	SET_STRING_ELT(classname, 0, mkChar("data.frame"));
	classgets(ans, classname);

	// sort it by variable, ego and alter
	SEXP R_fcall1, ordering, R_fcall2, ansnew;
	PROTECT(R_fcall1 = lang4(install("order"), col1, col3, col4));
	PROTECT(ordering = eval(R_fcall1, R_GlobalEnv));
	// now sort the data frame using [.data.frame and ordering
	PROTECT(R_fcall2 = lang4(install("[.data.frame"),
			ans, ordering, R_MissingArg));
	PROTECT(ansnew = eval(R_fcall2, R_GlobalEnv));

	UNPROTECT(18);

	if (sort)
	{
		return ansnew;
	}
	else
	{
		return ans;
	}
}
/**
 * Create a data-frame-plus from a chain. (prints nicely with PrintValue)
 */
SEXP getChainDFPlus(const Chain& chain, bool sort)
{
	SEXP main;
	PROTECT(main = getChainDF(chain, sort));

	SEXP initial;
	const vector<MiniStep *> & rMiniSteps = chain.rInitialStateDifferences();
	PROTECT(initial = getDFFromVector(rMiniSteps, false));

	SEXP is;
	PROTECT(is = install("initialStateDifferences"));
	setAttrib(main, is, initial);

	SEXP end;
	PROTECT(end = getDFFromVector(chain.rEndStateDifferences(), false));

	SEXP es;
	PROTECT(es = install("endStateDifferences"));
	setAttrib(main, es, end);

	SEXP classname;
	PROTECT(classname = allocVector(STRSXP, 2));
	SET_STRING_ELT(classname, 0, mkChar("chains.data.frame"));
	SET_STRING_ELT(classname, 1, mkChar("data.frame"));
	classgets(main, classname);

	UNPROTECT(6);
	return main;
}



/** Create a list from a ministep. Easy to create, but prints untidily!
 *
 */
SEXP getMiniStepList(const MiniStep& miniStep, int period)
{
	SEXP MINISTEP;
	PROTECT(MINISTEP = allocVector(VECSXP, 13));
	// unused elements (9, 10) used to contain the change contributions.
	SET_VECTOR_ELT(MINISTEP, 3, ScalarInteger(miniStep.ego()));
	if (miniStep.networkMiniStep())
	{
		const NetworkChange& networkChange =
			dynamic_cast<const NetworkChange &>(miniStep);
		SET_VECTOR_ELT(MINISTEP, 0, mkString("Network"));
		SET_VECTOR_ELT(MINISTEP, 4, ScalarInteger(networkChange.alter()));
		SET_VECTOR_ELT(MINISTEP, 5, ScalarInteger(0));
	}
	else
	{
		SET_VECTOR_ELT(MINISTEP, 0, mkString("Behavior"));
		const BehaviorChange& behaviorChange =
			dynamic_cast<const BehaviorChange &>(miniStep);
		SET_VECTOR_ELT(MINISTEP, 4, ScalarInteger(0));
		SET_VECTOR_ELT(MINISTEP,
			5,
			ScalarInteger(behaviorChange.difference()));
	}
	SET_VECTOR_ELT(MINISTEP, 1, ScalarInteger(miniStep.variableId()));
	SET_VECTOR_ELT(MINISTEP, 11, ScalarLogical(miniStep.missing(period)));
	SET_VECTOR_ELT(MINISTEP, 12, ScalarLogical(miniStep.diagonal()));
	SET_VECTOR_ELT(MINISTEP, 2, mkString(miniStep.variableName().c_str()));
	SET_VECTOR_ELT(MINISTEP, 7,	ScalarReal(miniStep.logOptionSetProbability()));
	SET_VECTOR_ELT(MINISTEP, 8, ScalarReal(miniStep.logChoiceProbability()));
	SET_VECTOR_ELT(MINISTEP, 6, ScalarReal(miniStep.reciprocalRate()));

	UNPROTECT(1);
	return MINISTEP;
}

/**
 * Create a list from a chain. Easy to create, but prints untidily!
 */
SEXP getChainList(const Chain& chain)
{
	SEXP ans;

	PROTECT(ans = allocVector(VECSXP, chain.ministepCount() - 1));

	MiniStep *pMiniStep = chain.pFirst()->pNext();
	for (int i = 0; i < chain.ministepCount() - 1; i++)
	{
		SET_VECTOR_ELT(ans, i, getMiniStepList(*pMiniStep, chain.period()));
		pMiniStep = pMiniStep->pNext();
	}

	// Add mu, sigma as attributes
	SEXP mu, sigma2, finalReciprocalRate;
	PROTECT(mu = allocVector(REALSXP, 1));
	REAL(mu)[0] = chain.mu();
	SEXP muu;
	PROTECT(muu = install("mu"));
	setAttrib(ans, muu, mu);
	PROTECT(sigma2 = allocVector(REALSXP, 1));
	REAL(sigma2)[0] = chain.sigma2();
	SEXP sigma;
	PROTECT(sigma = install("sigma2"));
	setAttrib(ans, sigma, sigma2);
	PROTECT(finalReciprocalRate = allocVector(REALSXP, 1));
	REAL(finalReciprocalRate)[0] = chain.finalReciprocalRate();
	SEXP frr;
	PROTECT(frr = install("finalReciprocalRate"));
	setAttrib(ans, frr, finalReciprocalRate);
	// get the initial state ministeps
	SEXP initial;
	int numberInitial = chain.rInitialStateDifferences().size();
	PROTECT(initial = allocVector(VECSXP, numberInitial));

	for (int i = 0; i < numberInitial; i++)
	{
		const MiniStep * pMiniStep2 = (chain.rInitialStateDifferences())[i];
		SET_VECTOR_ELT(initial, i,
			getMiniStepList(*pMiniStep2, chain.period()));
	}
	SEXP init;

	PROTECT(init = install("initialStateDifferences"));
	setAttrib(ans, init, initial);

	// get the end state ministeps
	SEXP end;
	int numberEnd = chain.rEndStateDifferences().size();
	PROTECT(end = allocVector(VECSXP, numberEnd));

	for (int i = 0; i < numberEnd; i++)
	{
		const MiniStep * pMiniStep2 = (chain.rEndStateDifferences())[i];
		SET_VECTOR_ELT(end, i,
			getMiniStepList(*pMiniStep2, chain.period()));
	}
	SEXP en;

	PROTECT(en = install("endStateDifferences"));
	setAttrib(ans, en, end);
	UNPROTECT(11);
	return ans;
}

/**
 * Create a list of tie flip contributions or behavior change contributions for each ministep in a chain
 */
SEXP getChangeContributionsList(const Chain& chain, SEXP EFFECTSLIST)
{
	// get the column names from the names attribute
	SEXP cols;
	PROTECT(cols = install("names"));
	SEXP Names = getAttrib(VECTOR_ELT(EFFECTSLIST, 0), cols);

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

		getColNos(Names, &netTypeCol, &nameCol, &effectCol,
				&parmCol, &int1Col, &int2Col, &initValCol,
				&typeCol, &groupCol, &periodCol, &pointerCol,
				&rateTypeCol, &intptr1Col, &intptr2Col, &intptr3Col,
				&settingCol);

	MiniStep * pMiniStep = chain.pFirst()->pNext();

	SEXP CHANGECONTRIBUTIONS;
	PROTECT(CHANGECONTRIBUTIONS = allocVector(VECSXP, chain.ministepCount() - 1));
	for (int m = 0; m < chain.ministepCount() - 1; m++)
	{
		NetworkChange * pNetworkChange = dynamic_cast<NetworkChange *>(pMiniStep);
		BehaviorChange * pBehaviorChange = dynamic_cast<BehaviorChange *>(pMiniStep);
		SEXP MINISTEPCONTRIBUTIONS = 0;
		SEXP EFFECTS;
		SEXP NETTYPE;
		PROTECT(NETTYPE = allocVector(STRSXP, 1));
		SEXP netType;
		PROTECT(netType = install("networkType"));
		if (pNetworkChange || pBehaviorChange)
		{
			const char * netwName;
			if(pNetworkChange)
			{
				netwName = pNetworkChange->variableName().c_str();
				SET_STRING_ELT(NETTYPE, 0, mkChar("oneMode"));
			}
			else
			{
				netwName = pBehaviorChange->variableName().c_str();
				SET_STRING_ELT(NETTYPE, 0, mkChar("behavior"));
			}
			for (int ii = 0; ii < length(EFFECTSLIST); ii++)
			{
				const char * networkName = CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, ii),nameCol), 0));
				if (strcmp(netwName, networkName) == 0)
				{
					SEXP NETNAME;
					PROTECT(NETNAME = allocVector(STRSXP, 1));
					SEXP netName;
					PROTECT(netName = install("networkName"));
					SET_STRING_ELT(NETNAME, 0, mkChar(networkName));
					EFFECTS = VECTOR_ELT(EFFECTSLIST, ii);
					map<const EffectInfo *, vector<double> >* contributions;
					if(pNetworkChange)
					{
						contributions = pNetworkChange->changeContributions();
					}
					else
					{
						contributions = pBehaviorChange->changeContributions();
					}
					int choices = contributions->begin()->second.size();
					int numberOfEffects =  length(VECTOR_ELT(EFFECTS,0));
					int rateEffects = 0;
					for(int e = 0; e < numberOfEffects; e++)
					{
						const char * effectType = CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), e));
						if (!(strcmp(effectType, "eval") == 0 || strcmp(effectType, "endow") == 0 || strcmp(effectType, "creation") == 0))
						{
							rateEffects = rateEffects + 1;
						}
					}
					int length = numberOfEffects-rateEffects;
					PROTECT(MINISTEPCONTRIBUTIONS = allocMatrix(REALSXP,length, choices));
					double * rcontr;
					rcontr = REAL(MINISTEPCONTRIBUTIONS);
					SEXP EFFECTNAMES;
					PROTECT(EFFECTNAMES = allocVector(STRSXP,length));
					SEXP effectNames;
					PROTECT(effectNames = install("effectNames"));
					SEXP EFFECTTYPES;
					PROTECT(EFFECTTYPES = allocVector(STRSXP,length));
					SEXP effectTypes;
					PROTECT(effectTypes = install("effectTypes"));
					int rates = 0;
					for (int i = 0; i < numberOfEffects; i++)
					{
						const char * effectType = CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
						if (strcmp(effectType, "eval") == 0 || strcmp(effectType, "endow") == 0 || strcmp(effectType, "creation") == 0)
						{
							EffectInfo * pEffectInfo = (EffectInfo *)R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(EFFECTS, pointerCol), i));
							SET_STRING_ELT(EFFECTNAMES, i-rates, mkChar(pEffectInfo->effectName().c_str()));
							SET_STRING_ELT(EFFECTTYPES, i-rates, mkChar(effectType));
							vector<double> values = (*contributions)[pEffectInfo];
							for(int a = 0; a < choices; a++)
							{
								rcontr[i - rates + a * length] = values.at(a);
							}
						}
						else
						{
							rates = rates+1;
						}
					}
					setAttrib(MINISTEPCONTRIBUTIONS, effectNames, EFFECTNAMES);
					setAttrib(MINISTEPCONTRIBUTIONS, effectTypes, EFFECTTYPES);
					setAttrib(MINISTEPCONTRIBUTIONS, netName, NETNAME);
					UNPROTECT(7);
				}
			}
		}
		setAttrib(MINISTEPCONTRIBUTIONS, netType, NETTYPE);
		SET_VECTOR_ELT(CHANGECONTRIBUTIONS, m, MINISTEPCONTRIBUTIONS);
		pMiniStep = pMiniStep->pNext();
		UNPROTECT(2);
	}
	UNPROTECT(2);
	return CHANGECONTRIBUTIONS;
}


SEXP createRObjectAttributes(SEXP EFFECTSLIST, SEXP& stats)
{
	int nEffects = 0;
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		nEffects += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}
	// get the column names from the names attribute
	SEXP cols;
	PROTECT(cols = install("names"));
	SEXP Names = getAttrib(VECTOR_ELT(EFFECTSLIST, 0), cols);

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

	int rateEffects = 0;
	vector<string> effNames;
	vector<string> effTypes;
	vector<string> netNames;
	vector<string> netTypes;

	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		for(int e = 0; e < length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0)); e++)
		{
			const char * effectType = CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), typeCol), e));
			if (strcmp(effectType, "eval") == 0 || strcmp(effectType, "endow") == 0 || strcmp(effectType, "creation") == 0)
			{
				EffectInfo * pEffectInfo = (EffectInfo *)R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), pointerCol), e));
				effNames.push_back(pEffectInfo->effectName());
				effTypes.push_back(effectType);
				netNames.push_back(CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), nameCol), e)));
				netTypes.push_back(CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), netTypeCol), e)));
			}
			else
			{
				rateEffects = rateEffects + 1;
			}
		}
	}
	int objEffects = nEffects-rateEffects;

	SEXP EFFECTNAMES = PROTECT(allocVector(STRSXP,objEffects));
	SEXP effectNames = PROTECT(install("effectNames"));
	SEXP EFFECTTYPES = PROTECT(allocVector(STRSXP,objEffects));
	SEXP effectTypes = PROTECT(install("effectTypes"));
	SEXP NETWORKNAMES = PROTECT(allocVector(STRSXP,objEffects));
	SEXP networkNames = PROTECT(install("networkNames"));
	SEXP NETWORKTYPES = PROTECT(allocVector(STRSXP,objEffects));
	SEXP networkTypes = PROTECT(install("networkTypes"));

	for(int eff = 0; eff < objEffects; eff++)
	{
		SET_STRING_ELT(EFFECTNAMES, eff, mkChar(effNames.at(eff).c_str()));
		SET_STRING_ELT(EFFECTTYPES, eff, mkChar(effTypes.at(eff).c_str()));
		SET_STRING_ELT(NETWORKNAMES, eff, mkChar(netNames.at(eff).c_str()));
		SET_STRING_ELT(NETWORKTYPES, eff, mkChar(netTypes.at(eff).c_str()));
	}
	if(stats)
	{
		setAttrib(stats, effectNames, EFFECTNAMES);
		setAttrib(stats, effectTypes, EFFECTTYPES);
		setAttrib(stats, networkNames, NETWORKNAMES);
		setAttrib(stats, networkTypes, NETWORKTYPES);
	}
	UNPROTECT(9);
	return NETWORKTYPES;
}

/** Create a ministep from a single ministep stored as a list (not dataframe).
 *
 */
MiniStep * makeMiniStepFromList(Data * pData, SEXP MINISTEP)
{
	if (strcmp(CHAR(STRING_ELT(VECTOR_ELT(MINISTEP, 0), 0)),
			"Network") == 0)
	{
		NetworkChange * pNetworkChange = new NetworkChange
			(pData->pNetworkData(CHAR(STRING_ELT(VECTOR_ELT(MINISTEP,
							2), 0))),
				asInteger(VECTOR_ELT(MINISTEP, 3)),
				asInteger(VECTOR_ELT(MINISTEP, 4)),
				asInteger(VECTOR_ELT(MINISTEP, 12)));
		return pNetworkChange;
	}
	else
	{
		BehaviorChange * pBehaviorChange = new BehaviorChange
			(pData->pBehaviorData(CHAR(STRING_ELT(VECTOR_ELT(MINISTEP,
							2), 0))),
				asInteger(VECTOR_ELT(MINISTEP, 3)),
				asInteger(VECTOR_ELT(MINISTEP, 5)));
		return pBehaviorChange;
	}
}

/**
 * Create a chain from a single chain stored as a list (not dataframe).
 */
Chain * makeChainFromList(Data * pData, SEXP CHAIN, int period)
{
	/* create a chain */
	Chain * pChain = new Chain(pData);

	/* set period */
	pChain->period(period);

	for (int i = 0; i < length(CHAIN); i++)
	{
		SEXP MINISTEP;
		MINISTEP = VECTOR_ELT(CHAIN, i);
		pChain->insertBefore(makeMiniStepFromList(pData, MINISTEP),
			pChain->pLast());
	}

    SEXP init;
    PROTECT(init = install("initialStateDifferences"));
    SEXP initialState = getAttrib(CHAIN, init);
	for (int i = 0; i < length(initialState); i++)
	{
		SEXP MINISTEP;
		MINISTEP = VECTOR_ELT(initialState, i);
		pChain->addInitialStateDifference(makeMiniStepFromList(pData,
				MINISTEP));
	}
	UNPROTECT(1);
	return pChain;
}

} // namespace siena
