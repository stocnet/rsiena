/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoOutDegreeFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * EgoOutDegreeFunction.
 *****************************************************************************/

#include "EgoOutDegreeFunction.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
EgoOutDegreeFunction::EgoOutDegreeFunction(string networkName) :
	NetworkAlterFunction(networkName)
{
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double EgoOutDegreeFunction::value(int alter)
{
	return this->pNetwork()->outDegree(this->ego());
}


/**
 * Returns the value of this function as an integer.
 */
int EgoOutDegreeFunction::intValue(int alter)
{
	return this->pNetwork()->outDegree(this->ego());
}

}
