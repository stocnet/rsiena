/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndirectTiesFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * IndirectTiesFunction.
 *****************************************************************************/
#include <stdexcept>

#include "IndirectTiesFunction.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/NetworkCache.h"
#include "utils/SqrtTable.h"



using namespace std;

namespace siena
{

/**
 * Constructor.
 */
IndirectTiesFunction::
	IndirectTiesFunction(std::string networkName,
				double parameter, bool firstIn, bool secondIn) :
	NetworkAlterFunction(networkName)
{
	this->lfirstin = firstIn;
	this->lsecondin = secondIn;
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = (parameter >= 2);
	this->lNbrDist2Nodes = 0;
}

/**
 * Deallocates this function.
 */
IndirectTiesFunction::~IndirectTiesFunction()
{
	delete[] this->lmark;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void IndirectTiesFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();
	int m = pNetwork->m();
	if (m > n)
	{
		n = m;
	}
	// now n will be OK for any choice of firstIn and secondIn.

	// A helper array of marks
	this->lmark = new bool[n];
	for (int i = 0; i < n; i++)
	{
		this->lmark[i] = false;
	}
}

/**
 * Does the necessary preprocessing work for calculating the
 * predicate for a specific ego. This method must be invoked before
 * calling IndirectTiesFunction::value(...).
 */

void IndirectTiesFunction::preprocessEgo(int ego)
{
	NetworkAlterFunction::preprocessEgo(ego);
	int statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();
	int m = pNetwork->m();
	if (m > n)
	{
		n = m;
	}
	IncidentTieIterator iterI;
	
	for (int i = 0; i < n; i++)
	{
		this->lmark[i] = false;
	}

	if (lfirstin)
	{
		iterI = pNetwork->inTies(ego);
	}
	else
	{
		iterI = pNetwork->outTies(ego);
	}

	// Invariant: mark[h] = true if and only if a two-path from ego
	// to h has been found.
	// Traverse all two-paths from ego

	for ( ;	iterI.valid(); iterI.next())
	{
		int k = iterI.actor();
		IncidentTieIterator iterJ;
		if (lsecondin)
		{
			iterJ = pNetwork->inTies(k);
		}
		else
		{
			iterJ = pNetwork->outTies(k);
		}

		for ( ; iterJ.valid(); iterJ.next())
		{
			int h = iterJ.actor();
			if ((!this->lmark[h]) && (h != ego))
			{
				// The first two-path from ego to h is found.
				this->lmark[h] = true;
				statistic ++;
			}
		}
	}
	
	if (pNetwork->isOneMode())
	{
// the direct ties have to be subtracted
		for ( ;	iterI.valid(); iterI.next())
		{
			if (this->lmark[iterI.actor()])
			{
				this->lmark[iterI.actor()] = false;
				statistic--;
			}			
		}		
	}
	
	this->lNbrDist2Nodes = statistic;
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */

double IndirectTiesFunction::value(int alter) const
{
	double value = this->lNbrDist2Nodes;
	if (this->lroot)
	{
		value = this->lsqrtTable->sqrt(value);
	}
	return value;
}

}
