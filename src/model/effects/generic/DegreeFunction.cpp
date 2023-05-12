/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DegreeFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * DegreeFunction.
 *****************************************************************************/

#include "DegreeFunction.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
DegreeFunction::DegreeFunction(string networkName, double par) :
	NetworkAlterFunction(networkName)
{
	this->lp = par;
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double DegreeFunction::value(int alter) const
{
	double statistic = 0;
	for (int j = 0; j < this->pNetwork()->n(); j++)
	{
		statistic = statistic + this->pNetwork()->outDegree(j);
	}
	return (statistic/this->pNetwork()->n()) - this->lp;
}

}
