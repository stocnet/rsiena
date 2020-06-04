/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HomCovariateMixedTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * HomCovariateMixedTwoPathFunction.
 *****************************************************************************/

#include <cmath>
#include "HomCovariateMixedTwoPathFunction.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "CovariateMixedNetworkAlterFunction.h"
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
HomCovariateMixedTwoPathFunction::HomCovariateMixedTwoPathFunction(
	string firstNetworkName, string secondNetworkName,
	string covariateName, bool excludeMissing) :
	CovariateMixedNetworkAlterFunction(firstNetworkName,
		secondNetworkName, covariateName)
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
void HomCovariateMixedTwoPathFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateMixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double HomCovariateMixedTwoPathFunction::value(int alter)
{
	int statistic = 0;
	if (!(this->lexcludeMissing && this->missing(alter)))
	{
		const Network * pFirstNetwork = this->pFirstNetwork();
		const Network * pSecondNetwork = this->pSecondNetwork();
		for (IncidentTieIterator iter = pSecondNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
			{
				// Get the receiver of the outgoing tie.
				int h = iter.actor();
				// 2-paths:
				if (!(this->lexcludeMissing && this->missing(h)))
					{
				if ((fabs(this->CovariateMixedNetworkAlterFunction::value(h) -
				this->CovariateMixedNetworkAlterFunction::value(this->ego()))
									< EPSILON) &&
				(fabs(this->CovariateMixedNetworkAlterFunction::value(alter)
				- this->CovariateMixedNetworkAlterFunction::value(this->ego()))
									< EPSILON) &&
								(pFirstNetwork->tieValue(h, alter) >= 1))
						{
							statistic++ ;
						}
					}
			}
	}
	return statistic;
}

}
