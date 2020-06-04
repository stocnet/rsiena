/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateFourCyclesEffect.cpp
 *
 * Description: This file contains the implementation of the
 * SameCovariateFourCyclesEffect class.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include "SameCovariateFourCyclesEffect.h"
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
SameCovariateFourCyclesEffect::SameCovariateFourCyclesEffect(const EffectInfo * pEffectInfo,
	bool TwoMode) :	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lTwoMode = TwoMode;
	this->lcounters = 0;

	if (pEffectInfo->internalEffectParameter() != 1 &&
		pEffectInfo->internalEffectParameter() != 2)
	{
		throw invalid_argument(
			"SameCovariateFourCyclesEffect: Parameter value 1 or 2 expected");
	}
	this->lroot = pEffectInfo->internalEffectParameter() == 2;
	this->lpSqrtTable = SqrtTable::instance();
}

/**
 * Destructor.
 */
SameCovariateFourCyclesEffect::~SameCovariateFourCyclesEffect()
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
void SameCovariateFourCyclesEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);

	delete[] this->lcounters;
	if (this->lTwoMode)
	{
		this->lcounters = new long int[this->pNetwork()->m()];
	}
	else
	{
		this->lcounters = new long int[this->pNetwork()->n()];
	}
}

/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void SameCovariateFourCyclesEffect::preprocessEgo(int ego)
{
	CovariateDependentNetworkEffect::preprocessEgo(ego);

	const Network * pNetwork = this->pNetwork();

	// Count the number of three paths i -> h <- k -> j from i to each j
	this->countThreePaths(ego, pNetwork, this->lcounters);

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
 * i -> h <- k -> j for v(k)=v(i)
 */
void SameCovariateFourCyclesEffect::countThreePaths(int i,
	const Network * pNetwork,
	long int * counters) const
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
	double myvalue = this->value(this->ego());

	// Enumerate all three-paths i -> h <- k -> j and update the counters.
	// The (average) time complexity is obviously O(d^3), where d is the
	// average degree.
	if (!(this->missing(this->ego())))
	{
		for (IncidentTieIterator iterI = pNetwork->outTies(i);
			iterI.valid();
			iterI.next())
		{
			int h = iterI.actor();

			for (IncidentTieIterator iterH = pNetwork->inTies(h);
				iterH.valid();
				iterH.next())
			{
				int k = iterH.actor();

				if ((i != k) && (!(this->missing(k))) &&
							(fabs(this->value(k) - myvalue) < EPSILON))
				{
					for (IncidentTieIterator iterK = pNetwork->outTies(k);
						iterK.valid();
						iterK.next())
					{
						int j = iterK.actor();

						if (j != h)
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
double SameCovariateFourCyclesEffect::calculateContribution(int alter) const
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
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double SameCovariateFourCyclesEffect::tieStatistic(int alter)
{
	// Avoid counting each 4-cycle four times in the evaluation statistic.
	// TODO: Is it okay to divide by 4 for endowment statistic as well?

	if (this->lroot)
	{
	return 0.5*(this->lpSqrtTable->sqrt(this->lcounters[alter]));
	}
	else
	{
	return this->lcounters[alter] * 0.25;
	}
}

}
