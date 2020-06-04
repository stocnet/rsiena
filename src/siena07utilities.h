/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07utilities.h
 *
 * Description: This file contains prototypes for various utilities used
 * for profiling, debugging, and creating R objects from C ones.
 *****************************************************************************/

#ifndef SIENA07UTILITIES_H_
#define SIENA07UTILITIES_H_

#include <Rinternals.h>

#include "network/Network.h"
#include "model/variables/DependentVariable.h"

namespace siena
{
	class Data;
	class BehaviorVariable;
	class Network;
	class MiniStep;
	class Chain;
	class State;
}

using namespace siena;

//--------------------------------------
// utility functions to process groups
//-------------------------------------

/* Calculate the period number of this group and period, to access
 * correct stored chain.
 */
int periodFromStart(std::vector<Data *> & pGroupData, int group, int period);

/* Calculate the total number of periods in all groups, which is the dimension
 * of some returned arrays.
 */

int totalPeriods(std::vector<Data *> & pGroupData);
/**
 * Traps errors so R can stop the function rather than being stoppped itself.
 *
 */
void Rterminate();

/**
 * print out the data for profiling with gprof
 *
 */
void printOutData(Data *pData);

SEXP getBehaviorValues(const BehaviorVariable & behavior);
SEXP getAdjacency(const Network& net);
SEXP getEdgeList(const Network& net);

/**
 * utilities to access chains and ministeps
 */
namespace siena
{
	SEXP var_to_sexp(DependentVariable * pVar);
	SEXP net_to_sexp(const Network * pVar);

	SEXP getMiniStepDF(const MiniStep & miniStep);
	SEXP getChainDF(const Chain& chain, bool sort=true);
	SEXP getChainDFPlus(const Chain & chain, bool sort=true);
	SEXP getDFFromVector(const std::vector<MiniStep *> & rMiniSteps, bool sort=true);
	SEXP getMiniStepList(const MiniStep & miniStep, int period);
	SEXP getChainList(const Chain & chain);
	SEXP getChangeContributionsList(const Chain & chain, SEXP EFFECTSLIST);
	SEXP createRObjectAttributes(SEXP EFFECTSLIST, SEXP & stats);
	Chain * makeChainFromList(Data * pData, SEXP CHAIN, int period);
	MiniStep * makeMiniStepFromList(Data * pData, SEXP MINISTEP);
	Chain * createMissingChain(int period, Data * data, const State& initialState);
}

#endif // SIENA07UTILITIES_H_
