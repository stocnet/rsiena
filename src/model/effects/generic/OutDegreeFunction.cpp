/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutDegreeFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * OutDegreeFunction.
 *****************************************************************************/

#include "OutDegreeFunction.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
OutDegreeFunction::OutDegreeFunction(string networkName) :
	NetworkAlterFunction(networkName)
{
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double OutDegreeFunction::value(int alter)
{
	return this->pNetwork()->outDegree(alter);
}


/**
 * Returns the value of this function as an integer.
 */
int OutDegreeFunction::intValue(int alter)
{
	return this->pNetwork()->outDegree(alter);
}

}
