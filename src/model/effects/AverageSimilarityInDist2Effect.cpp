/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageSimilarityInDist2Effect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageSimilarityInDist2Effect class.
 *****************************************************************************/

#include "AverageSimilarityInDist2Effect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "NetworkDependentBehaviorEffect.h"
#include "data/BehaviorLongitudinalData.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AverageSimilarityInDist2Effect::AverageSimilarityInDist2Effect(
	const EffectInfo * pEffectInfo, bool divide) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
}

/**
 * Deallocates this effect object;
 */
AverageSimilarityInDist2Effect::~AverageSimilarityInDist2Effect()
{
}

/**
 * Calculates the change in similarity with z_alt if adding 1 to z_ego,
 * except for division by range; this is done later.
 */
double AverageSimilarityInDist2Effect::changesim(double zalt, double zego) const
{
	double thechange = 0;
	if (zego > zalt)
	{
		thechange = -1;
	}
	else if (zego < zalt -1)
	{
		thechange = 1;
	}
	else
	{
		thechange = 2*(zalt - zego) - 1;
	}
	return thechange;
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 * It is assumed that preprocessEgo(ego) has been called before.
 */
double AverageSimilarityInDist2Effect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	if (!(difference == 0))
	{
		const Network * pNetwork = this->pNetwork();

		if (pNetwork->outDegree(actor) > 0)
		{
		// The formula for the effect:
		// s_i(x) = sum_{all neighbors j of i} {sim(v_j, z_i) - average sim},
		// where v_j is the average behavior of j's in-neighbors,
		// excluding i.
		// We need to calculate the change delta in s_i(x), if we changed
		// z_i to z_i + difference.
		// if (not divide), s_i(x) is further divided by outdegree(i).
			for (IncidentTieIterator iter = pNetwork->outTies(actor);
				iter.valid();
				iter.next())
			{
				int deg = pNetwork->inDegree(iter.actor());
				double alterValue = 0;
				if (deg > 0)
				{
					alterValue = this->totalInAlterValue(iter.actor());
// defined in NetworkDependentBehaviorEffect.cpp
					if (this->pNetwork()->tieValue(actor, iter.actor()) >= 1)
					{
						alterValue -= this->centeredValue(actor);
						deg--;
					}
				}
				if (deg > 0)
				{
					alterValue /= deg;
// now altervalue is $\check z_j^{(-i)}$ 
// defined in the manual for j = iter.actor()
					if (difference > 0) // then difference will be 1
					{
						contribution +=
							changesim(alterValue, this->centeredValue(actor));
					}
					else if (difference < 0)
					{
						contribution -=
						  changesim(alterValue, this->centeredValue(actor) - 1);
					}
				}
			}
			contribution /= this->range();
			if (this->ldivide)
			{
				contribution /= pNetwork->outDegree(actor);
			}
		}
	}
	return contribution;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * currentValues given for the behavior variable.
 */
double AverageSimilarityInDist2Effect::egoStatistic(int i, 
													double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(i) > 0)
	{
		for (IncidentTieIterator iter = pNetwork->outTies(i);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();
			int degree = this->pNetwork()->inDegree(j);
			if (degree > 0)
			{
				double inAlter = 0;
				for (IncidentTieIterator iterj = pNetwork->inTies(j);
				iterj.valid();
				iterj.next())
				{
					int k = iterj.actor();
					if (!(k == i))
					{
						inAlter += currentValues[k];
					}
					else
					{
						degree--;
					}
				}
				if (degree > 0)
				{
					inAlter /= degree;
					statistic += this->similarity(currentValues[i], inAlter);
				}
			}
		}
		if (this->ldivide)
		{
			statistic /= pNetwork->outDegree(i);
		}
	}
	return statistic;
}

}
