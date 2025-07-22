/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameInCovariateFourCyclesEffect.cpp
 *
 * Description: This file contains the implementation of the
 * SameInCovariateFourCyclesEffect class.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include "SameInCovariateFourCyclesEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
SameInCovariateFourCyclesEffect::SameInCovariateFourCyclesEffect(const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lTwoMode = false;
	this->lcounters = 0;

	if (pEffectInfo->internalEffectParameter() != 1 &&
		pEffectInfo->internalEffectParameter() != 2)
	{
		throw invalid_argument(
			"SameInCovariateFourCyclesEffect: Parameter value 1 or 2 expected");
	}

	this->lroot = pEffectInfo->internalEffectParameter() == 2;
	this->lpSqrtTable = SqrtTable::instance();
}


/**
 * Destructor.
 */
SameInCovariateFourCyclesEffect::~SameInCovariateFourCyclesEffect()
{
	delete[] this->lcounters;
	this->lcounters = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void SameInCovariateFourCyclesEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);
   this->lTwoMode = !this->pNetwork()->isOneMode();

	delete[] this->lcounters;
	
	int nm = this->pNetwork()->n();
	if (this->lTwoMode)
	{
		nm = this->pNetwork()->m();
	}

	this->lcounters = new long int[nm];
	for (int j = 0; j < nm; j++)
	{
		this->lcounters[j] = 0;
	}
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void SameInCovariateFourCyclesEffect::preprocessEgo(int ego)
{
	CovariateDependentNetworkEffect::preprocessEgo(ego);

	const Network * pNetwork = this->pNetwork();

 // Count the number of three paths i -> h <- k -> j for v(h)=v(j)
	this->countThreePaths(ego, pNetwork, this->lcounters, false);

	if (this->lroot)
	{
		// Count the number of 4-cycles the ego i is currently involved in.
		// This count is required for the sqrt case only.

		this->lcurrentCycleCount = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(ego);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();
			this->lcurrentCycleCount += this->lcounters[j];
		}

		// The above loop counted each 4-cycle twice
		this->lcurrentCycleCount /= 2;
	}
}


/**
 * For each j and the given i, this method calculates the number of three-paths
 * i -> h <- k -> j with v(h)=v(j)
 */
void SameInCovariateFourCyclesEffect::countThreePaths(int i,
	const Network * pNetwork, long int * counters, bool dropMissings) const
{
	// Initialize

	int nm = pNetwork->n();
	if (this->lTwoMode)
	{
		nm = pNetwork->m();
	}

	for (int j = 0; j < nm; j++)
	{
		counters[j] = 0;
	}

	// Enumerate all three-paths i -> h <- k -> j with v(h)=v(j) 
	// and update the counters.
	// The (average) time complexity is obviously O(d^3), where d is the
	// average degree.
	for (IncidentTieIterator iterI = pNetwork->outTies(i);
		iterI.valid();
		iterI.next())
	{
		int h = iterI.actor();
		
		if (!(this->missing(h) && dropMissings))
		{
			double hvalue = this->value(h);
			for (IncidentTieIterator iterH = pNetwork->inTies(h);
				iterH.valid();
				iterH.next())
			{
				int k = iterH.actor();

				if (i != k)
				{
					for (IncidentTieIterator iterK = pNetwork->outTies(k);
						iterK.valid();
						iterK.next())
					{
						int j = iterK.actor();
						if ((j != h) && (!(this->missing(j) && dropMissings)) &&
								((fabs(this->value(j) - hvalue) < EPSILON)))
						{
							counters[j]++;
						}
					}
				}
			}
		}
	}
}



/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double SameInCovariateFourCyclesEffect::calculateContribution(int alter) const
{
	double change;

	if (this->lroot)
	{
		int newCycleCount = this->lcurrentCycleCount;

		if (this->outTieExists(alter))
		{
			newCycleCount -= this->lcounters[alter];
		}
		else
		{
			newCycleCount += this->lcounters[alter];
		}

		change = fabs(this->lpSqrtTable->sqrt(newCycleCount) -
			this->lpSqrtTable->sqrt(this->lcurrentCycleCount));
	}
	else
	{
		change = this->lcounters[alter];
	}

	return change;
}


/**
 * The contribution of ego to the statistic.
 * Although preprocessEgo(ego) has been called before, it is not used,
 * because missing data have to be dealt with.
 * It would be more efficient to avoid 
 * counting each 4-cycle four times in the evaluation statistic.
 */
double SameInCovariateFourCyclesEffect::egoStatistic(int ego, const Network * pNetwork)
{
	double statistic = 0;
	
	if (pNetwork->outDegree(ego) > 0)
	{
 // Count the number of three paths ego -> h <- k -> j with v(h)=v(j)
		this->countThreePaths(ego, pNetwork, this->lcounters, true);
		for (IncidentTieIterator iterI = pNetwork->outTies(ego);
			iterI.valid();
			iterI.next())
		{
			if (this->lroot)
			{
				statistic += 0.5*(this->lpSqrtTable->sqrt(this->lcounters[iterI.actor()]));
			}
			else
			{
				statistic += 0.25*this->lcounters[iterI.actor()];
			}
		}
	}
	return statistic;
}


}
