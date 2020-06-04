/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * SameCovariateTwoPathFunction.
 *****************************************************************************/

#include <cmath>
#include "SameCovariateTwoPathFunction.h"
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

SameCovariateTwoPathFunction::SameCovariateTwoPathFunction(
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
void SameCovariateTwoPathFunction::initialize(const Data * pData,
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
double SameCovariateTwoPathFunction::value(int alter)
{
	int statistic = 0;
	if  (!(this->lexcludeMissing && this->missing(alter)))
	{
		const Network * pNetwork = this->pNetwork();
		// Iterate over outgoing ties in network W 
		for (IncidentTieIterator iter =
				pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
			{
				// Get the receiver of the outgoing tie.
				int h = iter.actor();
				// 2-paths:
				if (!(this->lexcludeMissing && this->missing(h)))
					{
					if ((fabs(this->CovariateNetworkAlterFunction::value(h)
				- this->CovariateNetworkAlterFunction::value(this->ego()))
									< EPSILON) &&
					(pNetwork->tieValue(h, alter) >= 1))
						{
							statistic++ ;
						}
					}
			}
	}
	return statistic;
}

}
