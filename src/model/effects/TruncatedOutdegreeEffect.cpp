/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TruncatedOutdegreeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * TruncatedOutdegreeEffect class.
 *****************************************************************************/


#include <stdexcept>
#include "TruncatedOutdegreeEffect.h"
#include "network/Network.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
TruncatedOutdegreeEffect::TruncatedOutdegreeEffect(
	const EffectInfo * pEffectInfo, bool right, bool outIso) : NetworkEffect(pEffectInfo)
{
	this->lOutIso = outIso;
	this->lc = 1;
	this->lright = right;

	if (this->lOutIso)
	{
		this->lc = 1;
	}
	else
	{
		this->lc = int(pEffectInfo->internalEffectParameter() + 0.01);
	}
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
double TruncatedOutdegreeEffect::calculateContribution(int alter) const
{
	double change = 0;

	// Current out-degree
	int d =	this->pNetwork()->outDegree(this->ego());

	if (lright) // outTrunc or outIso
	{
		if (this->outTieExists(alter))
		{
		// After a tie withdrawal, the new out-degree would be d-1, and
		// the new effect value would have decreased by 1 if d <= this->lc

			if (d <= this->lc)
			{
				if (this->lOutIso)
				{
					change = -1;
				}
				else
				{
					change = 1;
				}
			}
		}
		else
		{
		// When introducing a new tie, the new out-degree would be d+1, and
		// the new effect value would have increased by 1 if d < this->lc
			if (d < this->lc)
			{
				if (this->lOutIso)
				{
					change = -1;
				}
				else
				{
					change = 1;
				}
			}
		}
	}
	else // More
	{
		if (this->outTieExists(alter))
		{
		// After a tie withdrawal, the new out-degree would be d-1, and
		// the new effect value would have decreased by 1 if d > this->lc
			if (d > this->lc)
			{
				change = 1;
			}
		}
		else
		{
		// When introducing a new tie, the new out-degree would be d+1, and
		// the new effect value would have increased by 1 if d >= this->lc
			if (d >= this->lc)
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
double TruncatedOutdegreeEffect::egoStatistic(int ego,
	const Network * pNetwork)
{
	int statistic =	this->pNetwork()->outDegree(this->ego());

	if (this->lOutIso)
	{
		if (statistic <= 0)
		{
			statistic = 1;
		}
		else
		{
			statistic = 0;
		}
	}
	else
	{
		if (this->lright)
		{
			if (statistic > this->lc)
			{
				statistic = this->lc;
			}
		}
		else
		{
			if (statistic > this->lc)
			{
				statistic = statistic - this->lc;
			}
			else
			{
				statistic = 0;
			}
		}
	}
	return statistic;
}

}
