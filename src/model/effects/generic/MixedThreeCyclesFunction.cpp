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
							string secondNetworkName, double parameter) :
				MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = (parameter == 2)||(parameter == 4);
	this->lcenter = (parameter >= 3);
	this->lpFirstInStarTable = 0;
	this->lvariableName = firstNetworkName;
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
	if (lcenter)
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
}

/**
 * For each j and the given i, this method calculates
 * the number of mixed three-paths  i (W)-> h (W)<- k (X)-> j
 * where W = FirstNetwork, X = SecondNetwork, i = this->ego(), j=alter
 *
 * To generalize this allowing other directions and network choices:
 * see OutActDistance2Function.cpp for an example.
 */
double MixedThreeCyclesFunction::value(int alter)
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
	return statistic;
}

}
