/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateInStarFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * SameCovariateInStarFunction.
 *****************************************************************************/

#include <cmath>
#include "SameCovariateInStarFunction.h"
#include "NetworkAlterFunction.h"
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
 * @param[in] excludeMissing whether to exclude missing values
 */

SameCovariateInStarFunction::SameCovariateInStarFunction(
	string networkName, string covariateName, bool excludeMissing) :
	CovariateNetworkAlterFunction(networkName, covariateName)
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
void SameCovariateInStarFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateNetworkAlterFunction::initialize(pData, pState, period, pCache);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double SameCovariateInStarFunction::value(int alter) const
{
	int statistic = 0;
	if  (!(this->lexcludeMissing && this->missing(alter)))
	{
		const Network * pNetwork = this->pNetwork();
		// Iterate over outgoing ties in network W
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
				iter.valid();
				iter.next())
		{
			// Get the sender of the outgoing tie.
			int h = iter.actor();
			// in-2-stars:
			if (!(this->lexcludeMissing && this->missing(h)))
			{
				if ((fabs(this->CovariateNetworkAlterFunction::covvalue(h)
								- this->CovariateNetworkAlterFunction::covvalue(this->ego()))
							< EPSILON) &&
						(pNetwork->tieValue(alter, h) >= 1))
				{
					statistic++ ;
				}
			}
		}
	}
	return statistic;
}

}
