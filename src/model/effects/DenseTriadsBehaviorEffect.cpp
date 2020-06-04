/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DenseTriadsBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DenseTriadsBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>
#include "DenseTriadsBehaviorEffect.h"
#include "model/EffectInfo.h"
#include "network/OneModeNetwork.h"
#include "network/CommonNeighborIterator.h"
#include "network/IncidentTieIterator.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
DenseTriadsBehaviorEffect::DenseTriadsBehaviorEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldensity = (int) pEffectInfo->internalEffectParameter();
	this->lmark = 0;
	this->lbaseMark = 0;

	if (this->ldensity != 5 && this->ldensity != 6)
	{
		throw invalid_argument("Parameter value 5 or 6 expected.");
	}
}


/**
 * Destructor.
 */
DenseTriadsBehaviorEffect::~DenseTriadsBehaviorEffect()
{
	delete[] this->lmark;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DenseTriadsBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkDependentBehaviorEffect::initialize(pData, pState, period, pCache);
	delete[] this->lmark;
	int n = this->pNetwork()->n();
	this->lmark = new int[n];
	this->lbaseMark = 0;

	for (int i = 0; i < n; i++)
	{
		this->lmark[i] = 0;
	}
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double DenseTriadsBehaviorEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->denseTriadCount(actor);
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double DenseTriadsBehaviorEffect::egoStatistic(int ego, double * currentValues)
{
	return currentValues[ego] * this->denseTriadCount(ego);
}


/**
 * Returns the number of dense triads that the given actor is involved in.
 */
int DenseTriadsBehaviorEffect::denseTriadCount(int i)
{
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(
			"One-mode network expected in DenseTriadsBehaviorEffect");
	}

	this->lbaseMark += 2;

	// Count for each actor h the number of ties between i and h
	// (0, 1, or 2). We represent these number as:
	// mark[h] = baseMark + 2 if there are mutual ties between i and h,
	// mark[h] = baseMark + 1 if only one of the mutual ties is present,
	// mark[h] <= baseMark otherwise.

	for (IncidentTieIterator iter = pNetwork->inTies(i);
		iter.valid();
		iter.next())
	{
		this->lmark[iter.actor()] = this->lbaseMark + 1;
	}

	for (IncidentTieIterator iter = pNetwork->outTies(i);
		iter.valid();
		iter.next())
	{
		if (this->lmark[iter.actor()] <= this->lbaseMark)
		{
			this->lmark[iter.actor()] = this->lbaseMark + 1;
		}
		else
		{
			this->lmark[iter.actor()]++;
		}
	}

	// Now count the number of dense triads

	int count = 0;

	if (this->ldensity == 6)
	{
		for (CommonNeighborIterator iterI = pNetwork->reciprocatedTies(i);
			iterI.valid();
			iterI.next())
		{
			int j = iterI.actor();

			for (CommonNeighborIterator iterJ = pNetwork->reciprocatedTies(j);
				iterJ.valid();
				iterJ.next())
			{
				int h = iterJ.actor();

				if (this->lmark[h] == this->lbaseMark + 2)
				{
					count++;
				}
			}
		}

		// Each triad was counted twice
		count /= 2;
	}
	else
	{
		for (CommonNeighborIterator iterI = pNetwork->reciprocatedTies(i);
			iterI.valid();
			iterI.next())
		{
			int j = iterI.actor();
			IncidentTieIterator outIter = pNetwork->outTies(j);
			IncidentTieIterator inIter = pNetwork->inTies(j);

			while (outIter.valid())
			{
				while (inIter.valid() && inIter.actor() < outIter.actor())
				{
					inIter.next();
				}

				int h = outIter.actor();
				bool reciprocatedJH = inIter.valid() && inIter.actor() == h;
				int tieCountBetweenIH = this->lmark[h] - this->lbaseMark;

				if (reciprocatedJH)
				{
					// One more tie needed.

					if (tieCountBetweenIH > 0)
					{
						// Got the required tie. Make sure that the triad
						// is counted only once. If there's a reciprocated tie
						// between i and h, then we'll hit this triad twice,
						// but we count it only when j is the smaller of the
						// two alters.

						if (tieCountBetweenIH != 2 || j < h)
						{
							count++;
						}
					}
				}
				else
				{
					// Two more ties needed.

					if (tieCountBetweenIH == 2)
					{
						count++;
					}
				}

				outIter.next();
			}
		}
	}

	return count;
}

}
