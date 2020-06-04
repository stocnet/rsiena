/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocalDegreeContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ReciprocalDegreeContinuousEffect class.
 *****************************************************************************/

#include <cmath>
#include <string>
#include <stdexcept>
#include "ReciprocalDegreeContinuousEffect.h"
#include "network/OneModeNetwork.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
ReciprocalDegreeContinuousEffect::ReciprocalDegreeContinuousEffect(
	const EffectInfo * pEffectInfo, bool recip) :
		NetworkDependentContinuousEffect(pEffectInfo)
{
	this->lrecip = recip;
	// Indicates whether the effect addresses reciprocated or 
	// non-reciprocated ties
}


/**
 * Returns the number of reciprocated (or non-reciprocated) ties of a certain actor,
 * and thus how much this effect contributes to the change in the continuous behavior.
 */
double ReciprocalDegreeContinuousEffect::calculateChangeContribution(int actor)
{
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreeContinuousEffect");
	}

	double contribution = 0;
	if (lrecip)
	{
		contribution = sqrt(pNetwork->reciprocalDegree(actor));
	}
	else
	{
		contribution = sqrt(pNetwork->outDegree(actor) - pNetwork->reciprocalDegree(actor));
	}

	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double ReciprocalDegreeContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	return this->calculateChangeContribution(ego) * currentValues[ego];
}

}
