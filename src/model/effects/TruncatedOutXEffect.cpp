/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TruncatedOutXEffect.cpp
 *
 * Description: This file contains the implementation of the
 * TruncatedOutXEffect class.
 *****************************************************************************/

#include <stdexcept>
#include "TruncatedOutXEffect.h"
#include "network/Network.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
TruncatedOutXEffect::TruncatedOutXEffect(const EffectInfo * pEffectInfo) :
							CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lc = int(pEffectInfo->internalEffectParameter() + 0.001);
	// C++ always rounds downward
	if (this->lc < 1)
	{
		throw invalid_argument(
			"Truncated/More OutdegreeEffect: Parameter value must be at least 1");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TruncatedOutXEffect::calculateContribution(int alter) const
{
	double change = 0;
	const Network * pNetwork = this->pNetwork();

	if (this->value(alter) > 0)
	{
		int outd = -this->lc;
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
		{
			// Get the receiver of the outgoing tie.
			int h = iter.actor();
			if (this->value(h) > 0)
			{
				outd++;
			}
		}
// now outd =  \sum_h (x_{ih} I{v_h > 0}) - p.
		if (this->outTieExists(alter))
		{
			if (outd >= 1)
			{
				change = 1;
			}
		}
		else
		{
			if (outd >= 0)
			{
				change = 1;
			}
		}
	}
	return change;
}

/**
 * Calculates the statistic corresponding to the given ego. The parameter
 * pNetwork is always the current network as there are no endowment effects
 * of this kind.
 * TS: well, the endowment effect is implemented. I'm not sure about this.
 */
double TruncatedOutXEffect::egoStatistic(int ego,
	const Network * pSummationTieNetwork)
{
	double statistic = 0.0;
	const Network * pNetwork = this->pNetwork();
	int outd = -this->lc;
	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		iter.valid();
		iter.next())
	{
		// Get the receiver of the outgoing tie.
		int h = iter.actor();
		if (this->value(h) > 0)
		{
			outd++;
		}
	}
	if (outd > 0)
	{
		statistic = outd;
	}
	return statistic;
}
}
