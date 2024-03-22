/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutActDoubleDistance2Function.cpp
 *
 * Description: This file contains the implementation of the class
 * OutActDoubleDistance2Function.
 * It was built on OutActDoubleDistance2Function.cpp
 * and it might have been more efficient to make a double iterator class
 * and use that as iter in OutActDoubleDistance2Function.cpp.
 * If directionalities of the firstNetwork ties should be reversed,
 * look at OutActDoubleDistance2Function.cpp and reinstate firstIn here.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "OutActDoubleDistance2Function.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "data/Data.h"
#include "utils/SqrtTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */

OutActDoubleDistance2Function::OutActDoubleDistance2Function(
	string firstNetworkName, string secondNetworkName, double parameter,
					bool secondIn, bool average) :
					MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot =  (parameter >= 2);
	this->lsecondin = secondIn;
	this->laverage = average;
	this->lvariableName = secondNetworkName;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void OutActDoubleDistance2Function::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
	NetworkLongitudinalData * pNetworkData =
			pData->pNetworkData(this->lvariableName);
	if (!pNetworkData)
	{
		throw logic_error(
			"Network data for " + this->lvariableName + " expected.");
	}
	if (this->lsecondin)
	{
		this->lavdegree = pNetworkData->averageInDegree();
	}
	else
	{
		this->lavdegree = pNetworkData->averageOutDegree();
	}
	if (this->lroot)
	{
		this->lavdegree = sqrt(this->lavdegree);
		// not with lsqrtTable because not an integer
	}
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double OutActDoubleDistance2Function::value(int alter) const
{
	double statistic = 0;
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();
	IncidentTieIterator iter;
	int (Network::*pDegreeFunction)(int) const;
	int deg = 0;

	if (lsecondin)
	{
		pDegreeFunction = &siena::Network::inDegree;
	}
	else
	{
		pDegreeFunction = &siena::Network::outDegree;
	}

	for (iter = pFirstNetwork->outTies(this->ego()); iter.valid(); iter.next())
	{
		for (IncidentTieIterator iter2(pFirstNetwork->inTies(iter.actor()));
			iter2.valid();
			iter2.next())
			{
				if (this->lroot)
				{
					statistic +=
			(this->lsqrtTable->sqrt((pSecondNetwork->*pDegreeFunction)(iter2.actor()))
										- this->lavdegree);
					deg++;
				}
				else
				{
			statistic += ((pSecondNetwork->*pDegreeFunction)(iter2.actor())
														- this->lavdegree);
					deg++;
				}
			}
	}
	if ((deg > 0) && (this->laverage))
	{
		statistic /= deg;
	}
	return statistic;
}

}
