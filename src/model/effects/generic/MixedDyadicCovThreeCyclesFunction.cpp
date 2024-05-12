/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedDyadicCovThreeCyclesFunction.cpp
 *
 * Description: This file contains the implementation of the
 * MixedDyadicCovThreeCyclesFunction class.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include "MixedDyadicCovThreeCyclesFunction.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/tables/NetworkCache.h"
#include "network/CommonNeighborIterator.h"
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
MixedDyadicCovThreeCyclesFunction::MixedDyadicCovThreeCyclesFunction(string firstNetworkName,
			string secondNetworkName, string dyadicCovariateName, double parameter,
			bool average) :
		DyadicCovariateMixedNetworkAlterFunction(firstNetworkName, 
					secondNetworkName, dyadicCovariateName)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = ((parameter >= 4) && (parameter <= 6));
	this->lFirstWeight = ((parameter == 1)||(parameter == 3)||(parameter == 4)||(parameter == 6));
	this->lSecondWeight = ((parameter == 2)||(parameter == 3)||(parameter == 5)||(parameter == 6));
	this->lvariableName = firstNetworkName;
	this->laverage= average;
}

/**
 * Deallocates this function.
 */
MixedDyadicCovThreeCyclesFunction::~MixedDyadicCovThreeCyclesFunction()
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
void MixedDyadicCovThreeCyclesFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	DyadicCovariateMixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
	if ((this->firstNet_n()) != (this->secondNet_n()))
	{
		throw logic_error(
			"Two networks combined with different actor sets.");
	}
	if ((this->firstNet_n()) != (this->dyCov_n()))
	{
		throw logic_error(
			"Networks combined with dyadic covariate of incorrect first dimension.");
	}
	if ((this->firstNet_m()) != (this->dyCov_m()))
	{
		throw logic_error(
			"Networks combined with dyadic covariate of incorrect second dimension.");
	}
	this->ln = this->pFirstNetwork()->n();
	// A helper array of ltimesFound
	this->ltimesFound = new double[this->ln];
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

void MixedDyadicCovThreeCyclesFunction::preprocessEgo(int ego)
{
	MixedNetworkAlterFunction::preprocessEgo(ego);
	this->lsumDegs = 1;
	
	if (this->laverage)
	{	
		const Network * pFirstNetwork = this->pFirstNetwork();	
		double degs = 0;			
		this->lsumDegs = 0;
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
			int h = iter.actor();
			for (IncidentTieIterator iter2(pFirstNetwork->inTies(h));
				iter2.valid();
				iter2.next())
			{
				if (iter2.actor() != this->ego())
				{
					double mult = 1;
					if (this->lFirstWeight) 
					{
						mult *= this->dyadicValue(ego, h);
					}
					if (this->lSecondWeight) 
					{
						mult *= this->dyadicValue(iter2.actor(), h);
					}				
					if (this->lroot)
					{
						this->ltimesFound[iter2.actor()] += mult;
					}
					else
					{
						degs += mult;				
					}
				}
			}
		}
		if (this->lroot)
		{	
			degs = 0;
			for (int k = 0; k < this->ln; k++)
			{
				degs += sqrt(this->ltimesFound[k]);
			}
		}
		this->lsumDegs = degs;
	}
}

/**
 * For each j and the given i, this method calculates
 * the number of mixed three-paths  i (W)-> h (W)<- k (X)-> j
 * where W = FirstNetwork, X = SecondNetwork, i = this->ego(), j=alter;
 * W-paths are weighted by the dyadic covariate.
 *
 * To generalize this allowing other directions and network choices:
 * see OutActDistance2Function.cpp for an example.
 */
double MixedDyadicCovThreeCyclesFunction::value(int alter) const
{
	double statistic = 0;
	const Network * pSecondNetwork = this->pSecondNetwork();

	for (IncidentTieIterator iter = pSecondNetwork->inTies(alter);
		iter.valid();
		iter.next())
	{
		int h = iter.actor();
		if (!(h == this->ego()))
		{
			double contribution = 0;
			for (CommonNeighborIterator iterIH = this->firstNetworkInStars(this->ego(), h);
					iterIH.valid();
					iterIH.next())
			{
				int k = iterIH.actor();
				double contrib_k = 1;
				if (this->lFirstWeight) 
				{
					contrib_k *= this->dyadicValue(this->ego(), k);
				};
				if (this->lSecondWeight) 
				{
					contrib_k *= this->dyadicValue(h, k);
				}
				contribution += contrib_k;				
			}			
			if (this->lroot)
			{
				statistic += sqrt(contribution);
			}
			else
			{
				statistic += contribution;
			}
		}
	}
	if ((this->laverage) && (this->lsumDegs != 0))
	{
		statistic /= this->lsumDegs;
	}	
	return statistic;
}

}
