/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OneModeNetworkAlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * OneModeNetworkAlterFunction.
 *****************************************************************************/

#include <stdexcept>
#include "OneModeNetworkAlterFunction.h"
#include "network/OneModeNetwork.h"

using namespace std;

namespace siena
{

/**
 * Creates a new function depending on a one-mode network with the
 * given name.
 */
OneModeNetworkAlterFunction::OneModeNetworkAlterFunction(string networkName) :
	NetworkAlterFunction(networkName)
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void OneModeNetworkAlterFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);

	if (dynamic_cast<const OneModeNetwork *>(this->pNetwork()) == 0)
	{
		throw logic_error("One-mode network expected.");
	}
}

}
