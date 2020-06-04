/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MaxAlterContinuousEffect.cpp
 *
 * Description: This file contains the implementation of the
 * MaxAlterContinuousEffect class.
 *****************************************************************************/

#include <cmath>
#include "MaxAlterContinuousEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/ContinuousVariable.h"
namespace siena
{

/**
 * Constructor.
 */
MaxAlterContinuousEffect::MaxAlterContinuousEffect(
	const EffectInfo * pEffectInfo, bool minim) :
		NetworkDependentContinuousEffect(pEffectInfo)
{
	this->lminim = minim;
	// Indicates whether it will be a minimum instead of a maximum
}


/**
 * Returns the max (or min) of a certain actor's alters, and thus how
 * much this effect contributes to the change in the continuous behavior.
 */
double MaxAlterContinuousEffect::calculateChangeContribution(int actor)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
			if (lminim)
			{
				contribution = 1e8;
			}
			else
			{
				contribution = -1e8;
			}
		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			if (lminim)
			{
				if (this->centeredValue(iter.actor()) < contribution)
				{
					contribution = this->centeredValue(iter.actor());
				}
			}
			else
			{
				if (this->centeredValue(iter.actor()) > contribution)
				{
					contribution = this->centeredValue(iter.actor());
				}
			}
		}
	}
	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double MaxAlterContinuousEffect::egoStatistic(int ego, double * currentValues)
{
	return this->calculateChangeContribution(ego) * currentValues[ego];
}

}
