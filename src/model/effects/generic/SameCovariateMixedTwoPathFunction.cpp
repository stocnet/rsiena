/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateMixedTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * SameCovariateMixedTwoPathFunction.
 *****************************************************************************/

#include <cmath>
#include "SameCovariateMixedTwoPathFunction.h"
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

SameCovariateMixedTwoPathFunction::SameCovariateMixedTwoPathFunction(
	string firstNetworkName, string secondNetworkName,
	string covariateName, bool same, bool excludeMissing) :
	CovariateMixedNetworkAlterFunction(firstNetworkName,
		secondNetworkName, covariateName)
{
	this->lsame = same;
	this->lexcludeMissing = excludeMissing;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void SameCovariateMixedTwoPathFunction::initialize(const Data * pData,
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
double SameCovariateMixedTwoPathFunction::value(int alter) const
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
			if (this->lsame)
			{
				if (!(this->lexcludeMissing && this->missing(h)))
				{
					if ((fabs(this->CovariateMixedNetworkAlterFunction::covvalue(h) -
								this->CovariateMixedNetworkAlterFunction::covvalue(this->ego()))
							< EPSILON) &&
						(pFirstNetwork->tieValue(h, alter) >= 1))
					{
						statistic++ ;
					}
				}
			}
			else
			{
				if (!(this->lexcludeMissing && this->missing(h)))
				{
					if ((fabs(this->CovariateMixedNetworkAlterFunction::covvalue(h) -
								this->CovariateMixedNetworkAlterFunction::covvalue(this->ego()))
							>= EPSILON) &&
						(pFirstNetwork->tieValue(h, alter) >= 1))
					{
						statistic++ ;
					}
				}
			}
		}
	}
	return statistic;
}

}
