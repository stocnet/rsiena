/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleOutActFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * DoubleOutActFunction.
 *****************************************************************************/

#include <stdexcept>
#include "DoubleOutActFunction.h"
#include "network/Network.h"
#include "model/tables/TwoNetworkCache.h"
#include "network/CommonNeighborIterator.h"
#include "network/IncidentTieIterator.h"
#include "utils/SqrtTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
DoubleOutActFunction::DoubleOutActFunction(string firstNetworkName,
		string secondNetworkName, double parameter, bool change) :
	MixedNetworkAlterFunction(firstNetworkName, secondNetworkName)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = (parameter >= 2);
	this->lchange =  change;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DoubleOutActFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkAlterFunction::initialize(pData, pState, period, pCache);
	if (this->pFirstNetwork()->m() != this->pSecondNetwork()->m())
	{
		throw logic_error("doubleOutAct cannot be used for two-mode networks with different second modes.");
	}
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double DoubleOutActFunction::value(int alter)
{
	double statistic = 0;

	if (this->secondOutTieExists(alter))
	{
		const Network * pFirstNetwork = this->pFirstNetwork();
		const Network * pSecondNetwork = this->pSecondNetwork();
		for (CommonNeighborIterator iter(pFirstNetwork->outTies(this->ego()), pSecondNetwork->outTies(this->ego()));
				// for (IncidentTieIterator iter = pFirstNetwork->outTies(this->ego());
				iter.valid();
				iter.next())
		{
			// if (this->secondOutTieExists(iter.actor()))
			// {
			statistic++;
			// }
		}

		if (lchange)
		{
			if (this->lroot)
			{
				if (this->firstOutTieExists(alter))
				{
					statistic = (statistic * this->lsqrtTable->sqrt(statistic))
						- ((statistic - 1) * this->lsqrtTable->sqrt(statistic - 1));
				}
				else
				{
					statistic = ((statistic+1) * this->lsqrtTable->sqrt(statistic+1))
						- (statistic * this->lsqrtTable->sqrt(statistic));
				}
			}
			else
			{
				if (this->firstOutTieExists(alter))
				{
					statistic = 2*statistic - 1;
				}
				else
				{
					statistic = 2*statistic + 1;
				}
			}
		}
		else if (this->lroot)
		{
			statistic = (this->lsqrtTable->sqrt(statistic));
		}
	}

	return statistic;
}

}
