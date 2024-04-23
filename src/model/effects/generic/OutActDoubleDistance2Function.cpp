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
 * Deallocates this function.
 */
OutActDoubleDistance2Function::~OutActDoubleDistance2Function()
{
	delete[] this->ltimesFound;
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
	this->ln = this->pFirstNetwork()->n();
	// A helper array of ltimesFound
	this->ltimesFound = new int[this->ln] {};
	for (int h = 0; h < this->ln; h++)
	{
		this->ltimesFound[h] = 0;
	}
}


/**
 * increase is +1 if not lroot, and else sqrt(h+1)-sqrt(h)
 * used in pre-processing.
 */
double OutActDoubleDistance2Function::increase(int h) const
{
	double delta = 1;
	if (this->lroot)
	{
		delta = this->lsqrtTable->sqrt(h+1) - this->lsqrtTable->sqrt(h);
	}
	return delta;
}

/**
 * Does the necessary preprocessing work for calculating the
 * value for a specific ego. This method must be invoked before
 * calling OutActDoubleDistance2Function::value(...).
 */

void OutActDoubleDistance2Function::preprocessEgo(int ego)
{
	MixedNetworkAlterFunction::preprocessEgo(ego);
	double statistic = 0;
	
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();
	int (Network::*pDegreeFunction)(int) const;
	IncidentTieIterator iter;
	double degs = 0;
	for (int h = 0; h < this->ln; h++)
	{
		this->ltimesFound[h] = 0;
	}
	if (lsecondin)
	{
		pDegreeFunction = &siena::Network::inDegree;
	}
	else
	{
		pDegreeFunction = &siena::Network::outDegree;
	}

	// ltimesFound[h] is the number of times
	// that h was found at out-in distance 2 from ego
	// Traverse all out-in two-paths from ego

	for (iter = pFirstNetwork->outTies(ego); iter.valid(); iter.next())
	{
		for (IncidentTieIterator iter2(pFirstNetwork->inTies(iter.actor()));
			iter2.valid();
			iter2.next())
			{
				int h = iter2.actor();
				if (h != ego)
				{
					statistic += (this->increase(this->ltimesFound[h])) *
						((pSecondNetwork->*pDegreeFunction)(h) - this->lavdegree);
					degs += (this->increase(this->ltimesFound[h]));
					this->ltimesFound[h] ++;
				}
			}
	}
						
	if ((degs > 0) && (this->laverage))
	{
		statistic /= degs;
	}
	this->loutInDist2Degree = statistic;
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialize before and pre-processed with
 * respect to a certain ego.
 */
double OutActDoubleDistance2Function::value(int alter) const
{
	return this->loutInDist2Degree;
}

}
