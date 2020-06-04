/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutActDistance2Function.cpp
 *
 * Description: This file contains the implementation of the class
 * OutActDistance2Function.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "OutActDistance2Function.h"
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

OutActDistance2Function::OutActDistance2Function(
	string firstNetworkName, string secondNetworkName, double parameter,
					bool firstIn, bool secondIn) :
					MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot =  (parameter >= 2);
	this->lfirstin = firstIn;
	this->lsecondin = secondIn;
	this->lvariableName = secondNetworkName;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void OutActDistance2Function::initialize(const Data * pData,
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
double OutActDistance2Function::value(int alter)
{
	double statistic = 0;
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();
	IncidentTieIterator iter;
	int (Network::*pDegreeFunction)(int) const;

	if (lfirstin)
	{
		iter = pFirstNetwork->inTies(this->ego());
	}
	else
	{
		iter = pFirstNetwork->outTies(this->ego());
	}
	
	if (lsecondin)
	{
		pDegreeFunction = &siena::Network::inDegree;
	}
	else
	{
		pDegreeFunction = &siena::Network::outDegree;
	}
	
	for (  ; iter.valid(); iter.next())
			{
				if (this->lroot)
				{
					statistic +=
			(this->lsqrtTable->sqrt((pSecondNetwork->*pDegreeFunction)(iter.actor()))
										- this->lavdegree);
				}
				else
				{
			statistic += ((pSecondNetwork->*pDegreeFunction)(iter.actor()) 
														- this->lavdegree);
		}
	}
	return statistic;
}

}
