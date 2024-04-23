/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DegreeDistance2Function.cpp
 *
 * Description: This file contains the implementation of the class
 * DegreeDistance2Function.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "DegreeDistance2Function.h"
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

DegreeDistance2Function::DegreeDistance2Function(
	string networkName, double parameter,
					bool firstIn, bool secondIn, bool average) :
					NetworkAlterFunction(networkName)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot =  (parameter >= 2);
	this->lfirstin = firstIn;
	this->lsecondin = secondIn;
	this->laverage = average;
	this->lvariableName = networkName;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DegreeDistance2Function::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
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
double DegreeDistance2Function::value(int alter) const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	IncidentTieIterator iter;
	int (Network::*pDegreeFunction)(int) const;
	int deg = 0;

	if (lfirstin)
	{
		iter = pNetwork->inTies(this->ego());
	}
	else
	{
		iter = pNetwork->outTies(this->ego());
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
			(this->lsqrtTable->sqrt((pNetwork->*pDegreeFunction)(iter.actor()))
										- this->lavdegree);
					deg++;
				}
				else
				{
			statistic += ((pNetwork->*pDegreeFunction)(iter.actor()) 
														- this->lavdegree);
					deg++;
				}
	}
	if ((deg > 0) && (this->laverage))
	{
		statistic /= deg;		
	}
	return statistic;
}

}
