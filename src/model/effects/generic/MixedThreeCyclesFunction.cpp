/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedThreeCyclesFunction.cpp
 *
 * Description: This file contains the implementation of the
 * MixedThreeCyclesFunction class.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include "MixedThreeCyclesFunction.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"
#include "data/NetworkLongitudinalData.h"
#include "data/Data.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
MixedThreeCyclesFunction::MixedThreeCyclesFunction(string firstNetworkName,
							string secondNetworkName, double parameter,
							bool average) :
				MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = (parameter == 2)||(parameter == 4);
	this->lcenter = (parameter >= 3);
	this->lpFirstInStarTable = 0;
	this->lvariableName = firstNetworkName;
	this->laverage= average;
	if ((this->laverage) && (this->lcenter))
	{
		throw logic_error("sharedTo_Av can only have parameters 1 or 2");
	}
}

/**
 * Deallocates this function.
 */
MixedThreeCyclesFunction::~MixedThreeCyclesFunction()
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
void MixedThreeCyclesFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpFirstInStarTable = this->pFirstNetworkCache()->pInStarTable();
	NetworkLongitudinalData * pNetworkData =
		pData->pNetworkData(this->lvariableName);
	if (!pNetworkData)
	{
		throw logic_error("Network data for " + this->lvariableName + " expected.");
	}
	if (this->lcenter)
	{
		this->lavInTwoStar =
			(pNetworkData->averageSquaredInDegree() - pNetworkData->averageInDegree())
			/ (pNetworkData->m() - 1);
		if (this->lroot)
		{
			this->lavInTwoStar = sqrt(this->lavInTwoStar);
		}
	}
	else
	{
		this->lavInTwoStar = 0;
	}
	// A helper array of ltimesFound
	this->ln = this->pFirstNetwork()->n();
	this->ltimesFound = new int[this->ln];
	for (int h = 0; h < this->ln; h++)
	{
		this->ltimesFound[h] = 0;
	}
}


/**
 * Does the necessary preprocessing work for calculating the
 * value for a specific ego. This method must be invoked before
 * calling OutActDoubleDistance2Function::value(...).
 */

void MixedThreeCyclesFunction::preprocessEgo(int ego)
{
	MixedNetworkAlterFunction::preprocessEgo(ego);
	this->lsumDegs = 1;
	
	if (this->laverage) // prepare the denominator
	{
		double degs = 0;
		const Network * pFirstNetwork = this->pFirstNetwork();
		if (this->lroot)
		{
			for (int k = 0; k < this->ln; k++)
			{
				this->ltimesFound[k] = 0;
			}
		}
	
		for (IncidentTieIterator iter(pFirstNetwork->outTies(ego));
			iter.valid();
			iter.next())
		{	
			for (IncidentTieIterator iter2(pFirstNetwork->inTies(iter.actor()));
				iter2.valid();
				iter2.next())
			{
				if (iter2.actor() != this->ego())
				{
					if (this->lroot)
					{
						this->ltimesFound[iter2.actor()] ++;
					}
					else
					{
						degs ++;						
					}
				}
			}
		}
		if (this->lroot)
		{	
			degs = 0;
			for (int k = 0; k < this->ln; k++)
			{
				degs += this->lsqrtTable->sqrt(this->ltimesFound[k]);
			}
		}
		this->lsumDegs = degs;
	}
}

/**
 * For each j and the given i, this method calculates
 * the number of mixed three-paths  i (W)-> h (W)<- k (X)-> j
 * where W = FirstNetwork, X = SecondNetwork, i = this->ego(), j=alter
 *
 * To generalize this allowing other directions and network choices:
 * see OutActDistance2Function.cpp for an example.
 */
double MixedThreeCyclesFunction::value(int alter) const
{
	double statistic = 0;
	const Network * pSecondNetwork = this->pSecondNetwork();

	for (IncidentTieIterator iter = pSecondNetwork->inTies(alter);
		iter.valid();
		iter.next())
	{
		if (iter.actor() != this->ego())
		{
			if (this->lroot)
			{
				statistic +=
		(this->lsqrtTable->sqrt(this->lpFirstInStarTable->get(iter.actor())) -
							this->lavInTwoStar);
			}
			else
			{
				statistic +=
		(this->lpFirstInStarTable->get(iter.actor()) - this->lavInTwoStar);
			}
		}
	}
	if ((this->laverage) && (this->lsumDegs > 0))
	{
		statistic /= this->lsumDegs;
	}
	return statistic;
}

}
