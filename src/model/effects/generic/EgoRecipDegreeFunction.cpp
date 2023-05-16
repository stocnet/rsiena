/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoRecipDegreeFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * EgoRecipDegreeFunction.
 *****************************************************************************/

#include "EgoRecipDegreeFunction.h"
#include "network/OneModeNetwork.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
EgoRecipDegreeFunction::EgoRecipDegreeFunction(string networkName) :
	OneModeNetworkAlterFunction(networkName)
{
}



/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void EgoRecipDegreeFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	OneModeNetworkAlterFunction::initialize(pData, pState, period, pCache);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double EgoRecipDegreeFunction::value(int alter) const
{// but does not depend on alter
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());	
	return pONetwork->reciprocalDegree(this->ego());
}


/**
 * Returns the value of this function as an integer.
 */
int EgoRecipDegreeFunction::intValue(int alter)
{
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());	
	return pONetwork->reciprocalDegree(this->ego());
}

}
