/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InStarsTimesDegreesFunction.cpp
 *
 * Description: This file contains the implementation of the
 * InStarsTimesDegreesFunction class.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include "InStarsTimesDegreesFunction.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/tables/NetworkCache.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"
#include "data/Data.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
InStarsTimesDegreesFunction::InStarsTimesDegreesFunction(string firstNetworkName,
							string secondNetworkName, double parameter) :
				MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = (fabs(parameter - 2) < EPSILON);
	this->linv = (fabs(parameter + 1) < EPSILON);
//	this->lvariableName = secondNetworkName;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void InStarsTimesDegreesFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
}

/**
 * For each j and the given i, this method calculates
 * the number of two-in-stars i (W)-> h (W)<- j
 * weighted by the Z-indegree of h
 * where W = FirstNetwork, Z = SecondNetwork, i = this->ego(), j=alter.
 *
 * To generalize this allowing other directions and network choices:
 * see OutActDistance2Function.cpp for an example.
 */
double InStarsTimesDegreesFunction::value(int alter) const
{
	double statistic = 0;
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();

	CommonNeighborIterator iter(pFirstNetwork->outTies(this->ego()),
								pFirstNetwork->outTies(alter));

	for( ; iter.valid(); iter.next())
	{
		if (this->lroot)
		{
			statistic +=
			(this->lsqrtTable->sqrt(pSecondNetwork->inDegree(iter.actor())));
		}
		else
		{
			if (this->linv)
			{
				double indeg = pSecondNetwork->inDegree(iter.actor()) + 1;
				statistic += (1/indeg);
			}
			else
			{
			statistic += (pSecondNetwork->inDegree(iter.actor()));
			}
		}
	}
	return statistic;
}

}
