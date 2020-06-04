/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InteractionCovariateEffect.cpp
 *
 * Description: This file contains the implementation of the
 * InteractionCovariateEffect class.
 *****************************************************************************/

#include <stdexcept>
#include "InteractionCovariateEffect.h"
#include "model/EffectInfo.h"
#include "model/effects/SimilarityEffect.h"
#include "model/effects/AverageAlterEffect.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
InteractionCovariateEffect::InteractionCovariateEffect(
	const EffectInfo * pEffectInfo,
	bool averageSimilarity,
	bool totalSimilarity,
	bool averageAlter,
	bool totalAlter) :
		CovariateDependentBehaviorEffect(pEffectInfo)
{
	this->lpInteractingEffectInfo =
		new EffectInfo(pEffectInfo->variableName(),
			"Internal effect for interaction covariate effect",
			"",
			pEffectInfo->parameter(),
			0,
			pEffectInfo->interactionName2(),
			"",
			"");

	if (averageSimilarity)// || pEffectInfo->internalEffectParameter() == 1)
	{
		this->lpInteractingEffect =
			new SimilarityEffect(this->lpInteractingEffectInfo,
				true,
				false,
				false);
	}
	else if (totalSimilarity)// || pEffectInfo->internalEffectParameter() == 2)
	{
		this->lpInteractingEffect =
			new SimilarityEffect(this->lpInteractingEffectInfo,
				false,
				false,
				false);
	}
	else if (averageAlter) // || pEffectInfo->internalEffectParameter() == 3)
	{
		this->lpInteractingEffect =
			new AverageAlterEffect(this->lpInteractingEffectInfo, true, false);
	}
	else if (totalAlter) // || pEffectInfo->internalEffectParameter() == 4)
	{
		this->lpInteractingEffect =
			new AverageAlterEffect(this->lpInteractingEffectInfo, false, false);
	}
	else
	{
		// throw invalid_argument(
		// 	"Internal parameter should be in the range [1,3]");
		throw logic_error("Invalid call to Interaction Covariate Effect");
	}
}


/**
 * Destructor.
 */
InteractionCovariateEffect::~InteractionCovariateEffect()
{
	delete this->lpInteractingEffect;
	delete this->lpInteractingEffectInfo;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void InteractionCovariateEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentBehaviorEffect::initialize(pData,
		pState,
		period,
		pCache);
	this->lpInteractingEffect->initialize(pData, pState, period, pCache);
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double InteractionCovariateEffect::calculateChangeContribution(int actor,
	int difference)
{
	return this->covariateValue(actor) *
		this->lpInteractingEffect->calculateChangeContribution(actor,
			difference);
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double InteractionCovariateEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;

	// if (!this->missingCovariate(ego, this->period()) &&
	// 	!this->missingCovariate(ego, this->period() + 1))
	if (!this->missingCovariateEitherEnd(ego, this->period()))
	{
		statistic = this->covariateValue(ego) *
			this->lpInteractingEffect->egoStatistic(ego, currentValues);
	}

	return statistic;
}
/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double InteractionCovariateEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if (!this->missingCovariateEitherEnd(ego, this->period()))
	{
		statistic = this->covariateValue(ego) *
			this->lpInteractingEffect->egoEndowmentStatistic(ego, difference,
				currentValues);
	}

	return statistic;
}

}
