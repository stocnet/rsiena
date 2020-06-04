/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleInPopEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * DoubleInPopEffect.
 * Also see IndegreePopularityEffect.
 *****************************************************************************/

#include "DoubleInPopEffect.h"
#include "network/Network.h"
#include "network/CommonNeighborIterator.h"
#include "utils/SqrtTable.h"
#include "model/EffectInfo.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */

DoubleInPopEffect::DoubleInPopEffect(const EffectInfo * pEffectInfo) :
					MixedNetworkEffect(pEffectInfo)
{
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = (pEffectInfo->internalEffectParameter() >= 2);
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DoubleInPopEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkEffect::initialize(pData, pState, period, pCache);
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double DoubleInPopEffect::calculateContribution(int alter) const
{
	double statistic = 0;

	if (this->secondOutTieExists(alter))
	{
		statistic = 1;
		const Network * pFirstNetwork = this->pFirstNetwork();
		const Network * pSecondNetwork = this->pSecondNetwork();
		for (CommonNeighborIterator iter(pFirstNetwork->inTies(alter),
										pSecondNetwork->inTies(alter));
				iter.valid(); iter.next())
			{
				if (iter.actor() != this->ego())
				{
					statistic++;
				}
			}

		if (this->lroot)
		{
			statistic = this->lsqrtTable->sqrt(statistic);
		}
	}

	return statistic;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double DoubleInPopEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (this->secondOutTieExists(alter))
	{
		const Network * pFirstNetwork = this->pFirstNetwork();
		const Network * pSecondNetwork = this->pSecondNetwork();
		for (CommonNeighborIterator iter(pFirstNetwork->inTies(alter),
										pSecondNetwork->inTies(alter));
				iter.valid(); iter.next())
			{
				statistic++;
			}
		if (this->lroot)
		{
			statistic = (this->lsqrtTable->sqrt(statistic));
		}
	}

	return statistic;
}

}
