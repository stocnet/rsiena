/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InTieFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * InTieFunction.
 *****************************************************************************/

#include "InTieFunction.h"
#include "model/tables/NetworkCache.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
InTieFunction::InTieFunction(string networkName) :
	OneModeNetworkAlterFunction(networkName)
{
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double InTieFunction::value(int alter)
{
	return this->pNetworkCache()->inTieValue(alter);
}

}
