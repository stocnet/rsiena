/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedOnlyTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * MixedOnlyTwoPathFunction.
 *****************************************************************************/

#include <stdexcept>
#include "MixedOnlyTwoPathFunction.h"
#include "network/Network.h"
#include "model/tables/TwoNetworkCache.h"
#include "network/CommonNeighborIterator.h"
#include "network/IncidentTieIterator.h"
#include "utils/SqrtTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
MixedOnlyTwoPathFunction::MixedOnlyTwoPathFunction(string firstNetworkName,
		string secondNetworkName, double parameter) :
	MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void MixedOnlyTwoPathFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double MixedOnlyTwoPathFunction::value(int alter)
{
		double rvalue = this->intValue(alter);
		return rvalue;
}


/**
 * Returns the value of this function as an integer.
 */
int MixedOnlyTwoPathFunction::intValue(int alter)
{
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();
	bool tieThirdToAlter = false;
	bool noOtherTie = true;
	for (IncidentTieIterator iter1 = pSecondNetwork->outTies(this->ego());
			(iter1.valid() && noOtherTie);
			iter1.next())
	{
		int k = iter1.actor();  // "Third"
		for (IncidentTieIterator iter2 = pFirstNetwork->outTies(k);
			iter2.valid();
			iter2.next())
			if (iter2.actor() == alter)
			{
				tieThirdToAlter = true;				
			}
			else
			{
				noOtherTie = false;				
			}
	}
	if (tieThirdToAlter && noOtherTie)
	{
		return 1;
	}
	else		
	{
		return 0;
	}	
}

}
