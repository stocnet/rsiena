/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07models.h
 *
 * Description: This file contains prototypes for the siena simulation
 * modelling functions called from R
 *****************************************************************************/

#ifndef SIENA07MODELS_H_
#define SIENA07MODELS_H_

#include <Rinternals.h>

extern "C"
{

/**
 * Does one forward simulation for all the data by period within group
 */
SEXP forwardModel(SEXP DERIV, SEXP DATAPTR, SEXP SEEDS,
	SEXP FROMFINITEDIFF, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP RANDOMSEED2, SEXP RETURNDEPS, SEXP NEEDSEEDS,
	SEXP USESTREAMS, SEXP ADDCHAINTOSTORE, SEXP RETURNCHAINS, SEXP RETURNLOGLIK,
	SEXP RETURNACTORSTATISTICS, SEXP RETURNCHANGECONTRIBUTIONS);

/**
 * Does some MH steps for a specified group and period.
 * For multiple periods, the loop is always constructed in R.
 */
SEXP mlPeriod(SEXP DERIV, SEXP DATAPTR,
	SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP GROUP, SEXP PERIOD,
	SEXP NRUNMH, SEXP ADDCHAINTOSTORE, 
	SEXP RETURNDATAFRAME, SEXP RETURNDEPS, 
	SEXP RETURNCHAINS, SEXP RETURNLOGLIK, SEXP ONLYLOGLIK);


SEXP clearStoredChains(SEXP MODELPTR, SEXP KEEP, SEXP GROUPPERIOD);

SEXP getChainProbabilities(SEXP DATAPTR, SEXP MODELPTR,
	SEXP GROUP, SEXP PERIOD, SEXP INDEX, SEXP EFFECTSLIST, SEXP THETA,
	SEXP GETSCORES);

} // extern "C"

#endif /*SIENA07MODELS_H_*/
