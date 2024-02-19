/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SdeSimulation.cpp
 *
 * Description: This file contains the implementation of the
 * SdeSimulation class.
 *****************************************************************************/

#include <cmath> 
#include <R_ext/Error.h>
#include <R_ext/Print.h>

#include "EpochSimulation.h"
#include <R_ext/Error.h>
#include <R_ext/Print.h>
#include <Rinternals.h>
#include "EffectInfo.h"
#include "SdeSimulation.h"
#include "model/Model.h"
#include "model/effects/Effect.h"
#include "model/effects/ContinuousEffect.h"
#include "model/variables/ContinuousVariable.h"
#include "utils/Random.h"
#include <Rinternals.h>

namespace siena
{
/**
 * Creates a SDE simulation model for continuous variables.
 * @param pSimulation the owner simulation model of this model
 */
SdeSimulation::SdeSimulation(EpochSimulation * pSimulation)
{
	this->lpSimulation = pSimulation;
	int nContinuous = pSimulation->rContinuousVariables().size();
	this->lG = 1;
	this->lbasicScaleScore = 0;
	
	// this->A->resize(nContinuous, nContinuous);
	// this->G->resize(nContinuous, nContinuous);
	
	if (nContinuous > 1)
	{
		Rf_error("More than one continuous dependent variable: not implemented");
	}
	// for (unsigned i = 0; i < pSimulation->rcontinuousVariables.size(); i++) {

	const vector<Effect *> & rEffects = 
		pSimulation->rContinuousVariables()[0]->pFunction()->rEffects();
		
	for (unsigned i = 0; i < rEffects.size(); i++)
	{		
		ContinuousEffect * pEffect = (ContinuousEffect *) rEffects[i];
		if (pEffect->pEffectInfo()->effectName() == "feedback")
		{
			this->lA = pEffect->parameter();
		}
		else if (pEffect->pEffectInfo()->effectName() == "wiener")
		{
			this->lG = pEffect->parameter();
		}
	}
}

/**
 * Deallocates this object.
 */
SdeSimulation::~SdeSimulation()
{
	this->lpSimulation = 0;

}

/**
 * Initializes the object for a particular period.
 */
void SdeSimulation::initialize(int period)
{
	this->lperiod = period;
	this->lbasicScale = 
		this->lpSimulation->pModel()->basicScaleParameter(period);
}

/**
 * Computes the Bergstrom coefficients and stores these in the Effect objects.
 */
void SdeSimulation::setBergstromCoefficients(double dt)
{
	const vector<Effect *> & rEffects = 
		this->lpSimulation->rContinuousVariables()[0]->pFunction()->rEffects();

	this->lAdt = exp(this->lbasicScale * this->lA * dt);
	this->lQdt = (this->lAdt * this->lAdt - 1) * 
					this->lG * this->lG  / (2 * this->lA);
	double bdt = (this->lAdt - 1) / this->lA;
	
	for (unsigned i = 0; i < rEffects.size(); i++)
	{		
		ContinuousEffect * pEffect = (ContinuousEffect *) rEffects[i];
		if (pEffect->pEffectInfo()->effectName() == "feedback")
		{
			pEffect->coefficient(this->lAdt);
			// Rprintf("Adt = %f\n", this->lAdt);
		} 
		else if (pEffect->pEffectInfo()->effectName() == "wiener")
		{ 	
			// do nothing, the wiener effect does not contribute to the 
			// deterministic part of the Bergstrom formula
		}
		else
		{
			pEffect->coefficient(bdt * pEffect->parameter());
			// Rprintf("lBdt[%d] = %f\n", i, this->lBdt[i]);
		}
	}
}

// NYNKE: for multiple continuous variables this function will have to be
// 		  adapted: we then need one function that computes an n-dimensional
//		  random term and will have to use its elements for specific variables
double SdeSimulation::randomComponent() const
{
	return nextNormal(0, sqrt(this->lQdt));
}

double SdeSimulation::basicScaleScore() const
{
	return this->lbasicScaleScore;
}

void SdeSimulation::basicScaleScore(double score)
{
	this->lbasicScaleScore = score;
}

}
