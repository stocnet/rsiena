/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IsolateNetEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IsolateNetEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "IsolateNetEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
IsolateNetEffect::IsolateNetEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double IsolateNetEffect::calculateContribution(int alter) const
{
	double change = 0;

	if (this->pNetwork()->inDegree(this->ego()) == 0)
	{
		int degree = this->pNetwork()->outDegree(this->ego());
		if ((degree == 0) || ((degree == 1)&&(this->outTieExists(alter))))
		{
		// In the second case the single outtie is going to be withdrawn, 
		// so an isolate is created
			change = -1;
		}
	}

	return change;
}


/**
 * The contribution of ego to the statistic.
 * It is assumed that preprocessEgo(ego) has been called before.
 */
double IsolateNetEffect::egoStatistic(int ego, const Network * pNetwork)
{
	double statistic = 0;
	
	if ((this->pNetwork()->inDegree(this->ego()) == 0)
			&& (this->pNetwork()->outDegree(this->ego()) == 0))
	{
		statistic = 1;
	}

	return statistic;
}

}
