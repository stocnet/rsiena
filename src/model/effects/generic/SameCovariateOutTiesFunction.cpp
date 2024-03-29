/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateOutTiesFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * SameCovariateOutTiesFunction.
 *****************************************************************************/

#include <cmath>
#include "SameCovariateOutTiesFunction.h"
#include "NetworkAlterFunction.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
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

SameCovariateOutTiesFunction::SameCovariateOutTiesFunction(
	string networkName, string covariateName, bool sameValue, bool excludeMissing) :
	CovariateNetworkAlterFunction(networkName, covariateName)
{
	this->lsameValue = sameValue;
	this->lexcludeMissing = excludeMissing;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void SameCovariateOutTiesFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateNetworkAlterFunction::initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the
 * predicate for a specific ego. This method must be invoked before
 * calling SameCovariateOutTiesFunction::value(...).
 */

void SameCovariateOutTiesFunction::preprocessEgo(int ego)
{
	CovariateNetworkAlterFunction::preprocessEgo(ego);
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double SameCovariateOutTiesFunction::value(int alter) const
{
	int statistic = 0;
	if  (!(this->lexcludeMissing && this->missing(this->ego())))
	{
		const Network * pNetwork = this->pNetwork();
		// Iterate over incoming ties of alter
		if (lsameValue)
		{
			for (IncidentTieIterator iter =	pNetwork->outTies(alter);
			iter.valid();
			iter.next())
			{
				// Get the sender of the incoming tie.
				int h = iter.actor();
				// ego needs to have the same covariate value as h:
				if (!(this->lexcludeMissing && this->missing(h)))
				{
					if ((fabs(this->CovariateNetworkAlterFunction::covvalue(h)
					- this->CovariateNetworkAlterFunction::covvalue(this->ego()))
									< EPSILON))
					{
						statistic++ ;
					}
				}
			}
		}
		else
		{
			for (IncidentTieIterator iter =	pNetwork->inTies(alter);
			iter.valid();
			iter.next())
			{
				// Get the sender of the incoming tie.
				int h = iter.actor();
				// ego needs to have a different covariate value than h:
				if (!(this->lexcludeMissing && this->missing(h)))
				{
					if ((fabs(this->CovariateNetworkAlterFunction::covvalue(h)
				- this->CovariateNetworkAlterFunction::covvalue(this->ego()))
									>= EPSILON))
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
