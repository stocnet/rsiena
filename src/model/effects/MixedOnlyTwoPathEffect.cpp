/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedOnlyTwoPathEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * MixedOnlyTwoPathEffect.
 *****************************************************************************/

#include "MixedOnlyTwoPathEffect.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/CommonNeighborIterator.h"
#include "network/IncidentTieIterator.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
MixedOnlyTwoPathEffect::MixedOnlyTwoPathEffect(const EffectInfo * pEffectInfo) :
					MixedNetworkEffect(pEffectInfo)
{
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void MixedOnlyTwoPathEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	MixedNetworkEffect::initialize(pData, pState, period, pCache);
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pSecondNetwork());
	if (!pONetwork)
	{
		throw runtime_error(
			"One-mode second network expected in MixedOnlyTwoPathEffect");
	}
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double MixedOnlyTwoPathEffect::calculateContribution(int alter) const
{
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();
	double contribution = 0.0;

	for (IncidentTieIterator iter1 = pSecondNetwork->outTies(this->ego());
			iter1.valid();
			iter1.next())
	{
		int k = iter1.actor();  // "Third"
		bool ThirdTieToAlter = false;
		bool noOtherCommon = true;  // no other common tie to alter
		if (this->firstOutTieExists(alter))
		{
			for (CommonNeighborIterator iter2(pFirstNetwork->outTies(this->ego()),
					pFirstNetwork->outTies(k));
				(iter2.valid() && noOtherCommon);
				iter2.next())
			{
				if (iter2.actor() == alter)
				{
					ThirdTieToAlter = true;
				}
				else
				{
					noOtherCommon = false;
				}
			}
		}
		else
		{
			for (IncidentTieIterator iter2 = pFirstNetwork->outTies(k);
				(iter2.valid() && noOtherCommon);
				iter2.next())
			{
				int h = iter2.actor();
				if (h == alter)
				{
					ThirdTieToAlter = true;
				}
				else
				{
					if (this->secondOutTieExists(h))
					{
						noOtherCommon = false;
					}
				}
			}
		}
		if (noOtherCommon && ThirdTieToAlter)
		{
			contribution++;
		}
	}
	return contribution;
}

/**
 * Calculates the statistic corresponding to the given ego. The parameter
 * pNetwork is always the current network as there are no endowment effects
 * of this kind.
 */
double MixedOnlyTwoPathEffect::egoStatistic(int ego,
	const Network * pNetwork)
{
	const Network * pFirstNetwork = this->pFirstNetwork();
	const Network * pSecondNetwork = this->pSecondNetwork();
	double statistic = 0;
	for (IncidentTieIterator iter1 = pSecondNetwork->outTies(ego);
			iter1.valid();
			iter1.next())
	{
		int k = iter1.actor();  // "Third"
		bool notYet = true;
		for (CommonNeighborIterator iter2(pFirstNetwork->outTies(ego),
				pFirstNetwork->outTies(k));
				(iter2.valid() && notYet); iter2.next())
		{
			statistic++;
			notYet = false;
		}
	}
	return statistic;
}
}
