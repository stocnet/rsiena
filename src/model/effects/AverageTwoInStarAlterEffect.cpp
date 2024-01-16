/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageTwoInStarAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageTwoInStarAlterEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageTwoInStarAlterEffect.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/CommonNeighborIterator.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "NetworkDependentBehaviorEffect.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AverageTwoInStarAlterEffect::AverageTwoInStarAlterEffect(
	const EffectInfo * pEffectInfo, bool divide1, bool divide2) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldivide1 = divide1;
	// Indicates whether there will be division by the outdegree of ego (not used)
	this->ldivide2 = divide2;
	// Indicates whether there will be division by the indegree of alter (not used)
}

/**
 * Deallocates this effect object;
 */
AverageTwoInStarAlterEffect::~AverageTwoInStarAlterEffect()
{
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 * It is assumed that preprocessEgo(ego) has been called before.

 */
double AverageTwoInStarAlterEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	// is the if clause really necessary? The IncidentTieIterator should be able to handle this
	if (pNetwork->outDegree(actor) > 0) 
	{
		// The formula for the effect:
		// s_i(x) = v_i * sum(v_j) over all two-in-star alters j of i,
		// weighted by the number of common two-in-stars divided by
		// the sum of behavior of each two-in-star alter
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// divide1 and divide2 are not used for now

		double sumAlterValue = 0;
		double denom = 0;
		for (int j = 0; j < this->n(); j++)
		{
			if (j != actor)
			{
				double instars = 0;
				for (IncidentTieIterator iterH = pNetwork->outTies(j);
							iterH.valid();
							iterH.next())
					{
						int h = iterH.actor();
						for (IncidentTieIterator iterI = pNetwork->inTies(h);
									iterI.valid();
									iterI.next())
							{
								int i = iterI.actor();
								if (i == actor) // should stop after this - use while?
								{
									instars ++;
								}
							}
					}
				if (instars > 0)
				{
					denom += this->value(j);
				}
				sumAlterValue += this->value(j) * instars;
			}
		}
		contribution = difference * sumAlterValue;
		if (denom != 0)
		{
			contribution /= denom; //what happens if denom == 0 but sumaAlterValue != 0 ?
		}
	}
	return contribution;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * currentValues given for the behavior variable.
 */
double AverageTwoInStarAlterEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	double denom = 0;
	for (int j = 0; j < this->n(); j++)
	{
		if (j != ego)
		{
			double instarCount = 0;
			for (IncidentTieIterator iterH = pNetwork->outTies(j);
						iterH.valid();
						iterH.next())
				{
					int h = iterH.actor();
					for (IncidentTieIterator iterI = pNetwork->inTies(h);
								iterI.valid();
								iterI.next())
						{
							int i = iterI.actor();
							if (i == ego) // should stop after this - use while?
							{
								instarCount ++;
							}
						}
				}
			if (instarCount > 0)
			{
				denom += currentValues[j] + this->overallCenterMean();
			}
			statistic += (currentValues[j] + this->overallCenterMean()) * instarCount;
		}
	}
	statistic *= (currentValues[ego] + this->overallCenterMean());
	if (denom != 0)
	{
		statistic /= denom; //what happens if denom == 0 but statistic != 0 ?
	}
	
	return statistic;

}

}
