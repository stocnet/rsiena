/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InverseOutdegreeEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * InverseOutdegreeEffect.
 *****************************************************************************/

#include <stdexcept>
#include "InverseOutdegreeEffect.h"
#include "network/Network.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
InverseOutdegreeEffect::InverseOutdegreeEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lc = pEffectInfo->internalEffectParameter();

	if (this->lc < 1)
	{
		throw invalid_argument(
			"InverseOutdegreeEffect: Parameter value must be at least 1");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double InverseOutdegreeEffect::calculateContribution(int alter) const
{
	double sum =
		this->pNetwork()->outDegree(this->ego()) + this->lc;

	if (this->outTieExists(alter))
	{
		// Tie withdrawal
		return -1.0 / ((sum - 1) * sum);
	}
	else
	{
		// Tie introduction
		return -1.0 / ((sum + 1) * sum);
	}
}


/**
 * Calculates the statistic corresponding to the given ego. The parameter
 * pNetwork is always the current network as there are no endowment effects
 * of this kind.
 */
double InverseOutdegreeEffect::egoStatistic(int ego,
	const Network * pNetwork)
{
	return 1.0 / (pNetwork->outDegree(ego) + this->lc);
}

}
