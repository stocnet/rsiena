/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CatCovariateDependentNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDependentNetworkEffect class.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include "CatCovariateDependentNetworkEffect.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "network/Network.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
CatCovariateDependentNetworkEffect::CatCovariateDependentNetworkEffect(
		const EffectInfo * pEffectInfo) :
		CovariateDependentNetworkEffect(pEffectInfo), //
		lSimulatedOffset(0) {
	this->lpCovariateNumbers = 0;
}

/**
 * Constructor.
 *
 * @param pEffectInfo The effect info.
 * in addition, for gmom:
 * @param simulatedState If `true` the value(), missing() and actor_similarity()
 *        functions uses the simulated state, if any or the value at the end
 *        of the period.
 */
CatCovariateDependentNetworkEffect::CatCovariateDependentNetworkEffect(
		const EffectInfo * pEffectInfo, bool simulatedState) :
		CovariateDependentNetworkEffect(pEffectInfo), //
		lSimulatedOffset(simulatedState ? 1 : 0) {
	this->lpCovariateNumbers = 0;
}

/**
 * Destructor.
 */
CatCovariateDependentNetworkEffect::~CatCovariateDependentNetworkEffect()
{
	delete[] this->lpCovariateNumbers;
	this->lpCovariateNumbers = 0;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the state of the dependent variables at the beginning of the period
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CatCovariateDependentNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);
	
	int covMin = int(round(this->covariateMinimum()));
	int covMax = int(round(this->covariateMaximum()));
	if (covMin < 0)
	{
		throw logic_error("CatCovariateDependentNetworkEffect: minimum of first covariate is negative");		
	}
	if (covMax > 20)
	{
		throw logic_error("CatCovariateDependentNetworkEffect: first covariate has a maximum which is too large");		
	}
	
	covMax++;
	this->lpCovariateNumbers = new int[covMax] {};
	for (int i = 0; i < this->covarN(); i++)
	{
		this->lpCovariateNumbers[this->covariateIntValue(i)]++;		
	}
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the state of the dependent variables at the beginning of the period
 * @param[in] pSimulatedState the current simulated state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CatCovariateDependentNetworkEffect::initialize(const Data* pData,
		State* pState, State* pSimulatedState, int period, Cache* pCache) {
	CovariateDependentNetworkEffect::initialize(pData, pState, pSimulatedState, period, pCache);
	int covMin = int(round(this->covariateMinimum()));
	int covMax = int(round(this->covariateMaximum()));
	if (covMin < 0)
	{
		throw logic_error("CatCovariateDependentNetworkEffect: minimum of first covariate is negative");		
	}
	if (covMax > 20)
	{
		throw logic_error("CatCovariateDependentNetworkEffect: first covariate has a maximum which is too large");		
	}
	
	covMax++;
	this->lpCovariateNumbers = new int[covMax] {};
	int m = this->pNetwork()->m();

	for (int i = 0; i < m; i++)
	{
		this->lpCovariateNumbers[this->covariateIntValue(i)]++;		
	}
}


/**
 * Returns the covariate value for the given actor, rounded to integer.
 * For behavior, this is the non-centered value.
 */
int CatCovariateDependentNetworkEffect::covariateIntValue(int i) const
{
	return int(round(this->value(i)));
}


/**
 * Returns the numbers of nodes with rounded value a for the covariate.
 */
int CatCovariateDependentNetworkEffect::numberCovariateTies(int a) const
{
	return lpCovariateNumbers[a];
}


}
