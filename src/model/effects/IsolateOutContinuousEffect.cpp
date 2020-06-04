/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File:IsolateOutContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IsolateOutContinuousEffect class.
 *****************************************************************************/

#include <cmath>
#include "IsolateOutContinuousEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/ContinuousVariable.h"
namespace siena
{

/**
 * Constructor.
 */
IsolateOutContinuousEffect::IsolateOutContinuousEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentContinuousEffect(pEffectInfo)
{
}


/**
 * Returns the outdegree of a certain actor, and thus how much this effect
 * contributes to the change in the continuous behavior. 
 */
double IsolateOutContinuousEffect::calculateChangeContribution(int actor)
{
	double contribution = 0;
	
	if (this->pNetwork()->outDegree(actor) == 0) 
	{
		contribution = 1;
	}
	
	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double IsolateOutContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	
	if (this->pNetwork()->outDegree(ego) == 0) 
	{
		statistic = currentValues[ego];
	}
	
	return statistic;
}

}
