/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TruncatedOutdegreeEffect2.cpp
 *
 * Description: This file contains the implementation of the
 * TruncatedOutdegreeEffect2 class.
 *****************************************************************************/


#include <stdexcept>
#include "TruncatedOutdegreeEffect2.h"
#include "network/Network.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
TruncatedOutdegreeEffect2::TruncatedOutdegreeEffect2(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lc = pEffectInfo->internalEffectParameter();

	if (this->lc < 1)
	{
		throw invalid_argument(
			"TruncatedOutdegreeEffect: Parameter value must be at least 1");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TruncatedOutdegreeEffect2::calculateContribution(int alter) const
{
	double change = 0;

	// Current out-degree
	int d =	this->pNetwork()->outDegree(this->ego());

	if (this->outTieExists(alter))
	{
		// After a tie withdrawal, the new out-degree would be d-1, and
		// the new effect value would have decreased by 1 if d <= this->lc

		if (d <= this->lc)
		{
			change = 1;
		}
	}
	else
	{
		// When introducing a new tie, the new out-degree would be d+1, and
		// the new effect value would have increased by 1 if d < this->lc
		if (d < this->lc)
		{
			change = 1;
		}
	}

	return change;
}

/**
 * Calculates the statistic corresponding to the given ego. The parameter
 * pNetwork is always the current network as there are no endowment effects
 * of this kind.
 */
double TruncatedOutdegreeEffect2::egoStatistic(int ego,
	const Network * pNetwork)
{
	// Current out-degree
	int d =	this->pNetwork()->outDegree(this->ego());

	if (d <= this->lc)
	{
		return d;
	}
	else
	{
		return this->lc;
	}
}

}
