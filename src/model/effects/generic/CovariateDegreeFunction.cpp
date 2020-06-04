/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDegreeFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateDegreeFunction.
 *****************************************************************************/

#include <cmath>
#include "CovariateDegreeFunction.h"
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

CovariateDegreeFunction::CovariateDegreeFunction(
	string networkName, string covariateName, bool excludeMissing,
						bool incoming, bool forEgo, bool sqrtVersion):
	CovariateNetworkAlterFunction(networkName, covariateName)
{
	this->lexcludeMissing = excludeMissing;
	this->lincoming = incoming;
	this->lforEgo = forEgo;
	this->lsqrtVersion = sqrtVersion;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateDegreeFunction::initialize(const Data * pData,
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
double CovariateDegreeFunction::value(int alter)
{
	double statistic = 0;
	if  (!(this->lexcludeMissing && this->missing(alter)))
	{
		const Network * pNetwork = this->pNetwork();
		IncidentTieIterator iter;
		if (lincoming)
		{
			if (lforEgo)
			{
				iter = pNetwork->inTies(this->ego());
			}
			else
			{
				iter = pNetwork->inTies(alter);
			}
		}
		else
		{
			if (lforEgo)
			{
				iter = pNetwork->outTies(this->ego());
			}
			else
			{
				iter = pNetwork->outTies(alter);
			}
		}
		// iterate over ties in network W
		while (iter.valid())
			{
				// Get the actor at the other end of the tie.
				int h = iter.actor();
				if (!(this->lexcludeMissing && this->missing(h)))
					{
		statistic = statistic + this->CovariateNetworkAlterFunction::value(h);
					}
			iter.next();
			}
		if (lsqrtVersion)
		{
			if (statistic < 0)
			{
throw logic_error("param. 2 for mixed degree effect: only for covariate >= 0.");
			}
			statistic = sqrt(statistic);
		}
	}
	return statistic;
}
}
