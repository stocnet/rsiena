/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: XXWClosureEffect.cpp
 *
 * Description: This file contains the implementation of the
 * XXWClosureEffect class.
 *****************************************************************************/

#include "XXWClosureEffect.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "network/IncidentTieIterator.h"
#include "data/DyadicCovariateValueIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
XXWClosureEffect::XXWClosureEffect(const EffectInfo * pEffectInfo, bool outst, bool inst) :
	DyadicCovariateDependentNetworkEffect(pEffectInfo)
{
	this->loutStarSums = 0;
	this->linStarSums = 0;
	this->loutst = outst;
	this->linst = inst;
	this->lpNetwork = 0;
}


/**
 * Destructor.
 */
XXWClosureEffect::~XXWClosureEffect()
{
	delete[] this->loutStarSums;
	delete[] this->linStarSums;
	this->loutStarSums = 0;
	this->linStarSums = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void XXWClosureEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	DyadicCovariateDependentNetworkEffect::initialize(pData,
		pState, period, pCache);
	delete[] this->loutStarSums;
	delete[] this->linStarSums;
	this->loutStarSums = new double[this->pNetwork()->n()]; // only one-mode
	this->linStarSums = new double[this->pNetwork()->n()];
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void XXWClosureEffect::preprocessEgo(int ego)
{
	DyadicCovariateDependentNetworkEffect::preprocessEgo(ego);
	this->calculateOutStarSums(ego, this->pNetwork(), this->loutStarSums);
	this->calculateInStarSums(ego, this->pNetwork(), this->linStarSums);
}


/**
 * For a fixed i, this variable stores the value of sum_h x_{hi} w_{hj} for
 * each j.
 */
void XXWClosureEffect::calculateOutStarSums(int i,
	const Network * pNetwork, double * sums) const
{
	int n = pNetwork->n();

	// Initialize

	for (int j = 0; j < n; j++)
	{
		sums[j] = 0;
	}

	// Iterate over all h with x_{ih} = 1

	for (IncidentTieIterator iterH = pNetwork->inTies(i);
		iterH.valid();
		iterH.next())
	{
		int h = iterH.actor();

		// Iterate over all j with non-zero non-missing w_{hj}

		for (DyadicCovariateValueIterator iterJ = this->rowValues(h);
			iterJ.valid();
			iterJ.next())
		{
			int j = iterJ.actor();

			// Add the term x_{hi} w_{hj} (= w_{hj})
			sums[j] += iterJ.value();
		}
	}
}

/**
 * For a fixed i, this variable stores the value of sum_h w_{ih} x_{jh}
 * for each j.
 * note that this function differs from XWXClosureEffect.calculateInStarSums.
 */
void XXWClosureEffect::calculateInStarSums(int i,
	const Network * pNetwork,
	double * sums) const
{
	int n = pNetwork->n();

	// Initialize

	for (int j = 0; j < n; j++)
	{
		sums[j] = 0;
	}

	// Iterate over all h with non-zero non-missing w_{ih}

	for (DyadicCovariateValueIterator iterH = this->rowValues(i);
		iterH.valid();
		iterH.next())
		{
			int h = iterH.actor();

			for (IncidentTieIterator iterJ = pNetwork->inTies(h);
			iterJ.valid();
			iterJ.next())
			// Add the term w_{ih} x_{jh} (= w_{ih})
			{
				int j = iterJ.actor();
				sums[j] += iterH.value();
			}
		}
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double XXWClosureEffect::calculateContribution(int alter) const
{
	double statistic = 0;
	if (this->loutst) statistic = this->loutStarSums[alter];
	if (this->linst) statistic += this->linStarSums[alter];
	return statistic;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double XXWClosureEffect::tieStatistic(int alter)
{
	return this->loutStarSums[alter];
}

}
