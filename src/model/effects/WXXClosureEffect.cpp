/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WXXClosureEffect.cpp
 *
 * Description: This file contains the implementation of the
 * WXXClosureEffect class.
 *****************************************************************************/

#include "WXXClosureEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "data/DyadicCovariateValueIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
WXXClosureEffect::WXXClosureEffect(const EffectInfo * pEffectInfo) :
	DyadicCovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsums = 0;
}


/**
 * Destructor.
 */
WXXClosureEffect::~WXXClosureEffect()
{
	delete[] this->lsums;
	this->lsums = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void WXXClosureEffect::initialize(const Data * pData,
	State * pState, int period, Cache * pCache)
{
	DyadicCovariateDependentNetworkEffect::initialize(pData,
		pState, period, pCache);

	delete[] this->lsums;
	this->lsums = new double[this->pNetwork()->n()];
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void WXXClosureEffect::preprocessEgo(int ego)
{
	DyadicCovariateDependentNetworkEffect::preprocessEgo(ego);
	this->calculateSums(ego, this->pNetwork(), this->lsums);
}


/**
 * For each j and the given i, this method calculates the sum
 * sum_h w_{ih} x_{hj}.
 */
void WXXClosureEffect::calculateSums(int i,
	const Network * pNetwork, double * sums) const
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

		// Iterate over all j with a tie from h

		for (IncidentTieIterator iterJ = pNetwork->outTies(h);
				iterJ.valid();
				iterJ.next())
		{
			int j = iterJ.actor();

			// Add the term w_{ih} x_{hj} (= w_{ih})
			sums[j] += iterH.value();
		}
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double WXXClosureEffect::calculateContribution(int alter) const
{
	return this->lsums[alter];
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double WXXClosureEffect::tieStatistic(int alter)
{
	return this->lsums[alter];
}

}
