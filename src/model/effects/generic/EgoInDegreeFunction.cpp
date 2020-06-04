/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoInDegreeFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * EgoInDegreeFunction.
 *****************************************************************************/

#include "EgoInDegreeFunction.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
EgoInDegreeFunction::EgoInDegreeFunction(string networkName) :
	OneModeNetworkAlterFunction(networkName)
{
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double EgoInDegreeFunction::value(int alter)
{
	return this->pNetwork()->inDegree(this->ego());
}


/**
 * Returns the value of this function as an integer.
 */
int EgoInDegreeFunction::intValue(int alter)
{
	return this->pNetwork()->inDegree(this->ego());
}

}
