/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreeContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutdegreeEffect class.
 *****************************************************************************/

#include "OutdegreeContinuousEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/ContinuousVariable.h"
namespace siena
{

/**
 * Constructor.
 */
OutdegreeContinuousEffect::OutdegreeContinuousEffect(
	const EffectInfo * pEffectInfo, bool root) :
		NetworkDependentContinuousEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Returns the outdegree of a certain actor, and thus how much this effect
 * contributes to the change in the continuous behavior. 
 */
double OutdegreeContinuousEffect::calculateChangeContribution(int actor)
{
double change = this->pNetwork()->outDegree(actor);

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(change);
	}

	return change;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the continuous behavior variable.
 */
double OutdegreeContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = this->pNetwork()->inDegree(ego);
	
	if (this->lroot)
	{
		statistic = this->lsqrtTable->sqrt(statistic);
	}

	return currentValues[ego] * statistic;
}

}
