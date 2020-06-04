/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IsolatePopEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IsolatePopEffect class.
 *****************************************************************************/


#include "IsolatePopEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
IsolatePopEffect::IsolatePopEffect(
	const EffectInfo * pEffectInfo, bool outgoing) : NetworkEffect(pEffectInfo)
{
	this->loutgoing = outgoing;
	// Indicates whether the outdegree=0 of alter also is taken into account
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double IsolatePopEffect::calculateContribution(int alter) const
{
	double change = 0;

	if (this->loutgoing)
	{
		if (this->pNetwork()->outDegree(alter) == 0)
		{
			int degree = this->pNetwork()->inDegree(alter);
			if ((degree == 0) || ((degree == 1)&&(this->outTieExists(alter))))
			{
				// In the second case the single tie to this alter is going to be 
				// withdrawn, so an isolate is created
				change = 1;
			}
		}
	}
	else  // cannot be combined with the above in one if statement,
			// because for !outgoing this effect is also intended for bipartite networks
			// and then outDegree(alter) runs into an error.
	{
		int degree = this->pNetwork()->inDegree(alter);
		if ((degree == 0) || ((degree == 1)&&(this->outTieExists(alter))))
		{
			change = 1;
		}
	}
	return change;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double IsolatePopEffect::tieStatistic(int alter)
{
	double statistic = 0;
	
	if (this->loutgoing)
	{
		if ((this->pNetwork()->outDegree(alter) == 0)
				&& (this->pNetwork()->inDegree(alter) == 1))
		{
			statistic = 1;
		}
	}
	else // see above
	{
		if (this->pNetwork()->inDegree(alter) == 1)
		{
			statistic = 1;
		}
	}
	return statistic;
}


}
