/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File:WeightedMixedTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 *WeightedMixedTwoPathFunction.
 *****************************************************************************/

#include "WeightedMixedTwoPathFunction.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "utils/Utils.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 * @param[in] covariateName the name of the covariate this function is
 * associated with
 * @param[in] excludeMissing: whether to exclude missing values
 */

WeightedMixedTwoPathFunction::WeightedMixedTwoPathFunction(
	string firstNetworkName, string secondNetworkName,
	string dyadicCovariateName, bool excludeMissing) :
	DyadicCovariateMixedNetworkAlterFunction(firstNetworkName,
		secondNetworkName, dyadicCovariateName)
{
	this->lexcludeMissing = excludeMissing;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void WeightedMixedTwoPathFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	DyadicCovariateMixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double WeightedMixedTwoPathFunction::value(int alter)
{
	double statistic = 0;
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();
	for (IncidentTieIterator iter = pSecondNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
	{
		// Get the receiver of the outgoing tie.
		int h = iter.actor();
		// 2-paths:
		if (!(this->lexcludeMissing && this->missing(this->ego(),h)))
		{
			if (pFirstNetwork->tieValue(h, alter) >= 1)
			{
				statistic = statistic + this->dyadicValue(this->ego(), h);
			}
		}
	}
	return statistic;
}

}
