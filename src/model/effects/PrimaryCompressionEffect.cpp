/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PrimaryCompressionEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * PrimaryCompressionEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include <R_ext/Print.h>
#include "NetworkEffect.h"
#include "PrimaryCompressionEffect.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
PrimaryCompressionEffect::PrimaryCompressionEffect(const EffectInfo * pEffectInfo,
		bool inside, bool useSize): NetworkWithPrimaryEffect(pEffectInfo)
{
	this->llogNonPrimary = 0.0;
	this->lparameter = pEffectInfo->internalEffectParameter();
	this->luseSize = useSize;
	this->linside = inside;
	this->llogNonPrimary = 0;
	this->llogPrimary = 0;
	if ((this->lparameter <= 0) && (this->luseSize))
	{
		throw runtime_error(
	"Internal effect parameter should be positive in PrimaryCompressionEffect");
	}
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void PrimaryCompressionEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkWithPrimaryEffect::initialize(pData, pState, period, pCache);
}

/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void PrimaryCompressionEffect::preprocessEgo(int ego)
{
	NetworkWithPrimaryEffect::preprocessEgo(ego);
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();

	if (this->primaryDegree() < n-1)
	{
		double theRatio = (n - 1 - this->primaryDegree())/this->lparameter;
		if (theRatio < 1)
		{
			theRatio = 1.0;
		}
		this->llogNonPrimary = log(theRatio);
	}
	else // nothing outside primary setting of ego
	{
		this->llogNonPrimary = 0;
	}

	if (pNetwork->outDegree(ego) < this->primaryDegree())
	{
		double theRatio = (this->primaryDegree() - pNetwork->outDegree(ego))/this->lparameter;
		if (theRatio < 1)
		{
			theRatio = 1.0;
		}
		this->llogPrimary = log(theRatio);
	}
	else // ego has ties to all in ego's primary setting
	{
		this->llogPrimary = 0;
	}
}

/**
 * Returns if this effect is an ego effect.
 */
bool PrimaryCompressionEffect::egoEffect() const
{
	return true;
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double PrimaryCompressionEffect::calculateContribution(int alter) const
{
	double contribution = 0;
	if (this->linside)
	{
		if (inPrimarySetting(alter))
		{
			if (this->luseSize)
			{
				contribution = - this->llogPrimary;
			}
			else
			{
				contribution = 1;
			}
		}
	}
	else
	{
		if (!inPrimarySetting(alter))
		{
			if (this->luseSize)
			{
				contribution = - this->llogNonPrimary;
			}
		}
	}
	return contribution;
}

/**
 * The contribution of ego to the statistic.
 * For useSize, this effect is meant to be used a an offset, i.e., with fixed parameter.
 * Therefore this is meant to be used only if inside && !useSize
 */
double PrimaryCompressionEffect::egoStatistic(int ego, const Network * pNetwork)
{
	double statistic = 0;
//	int pd = this->primaryDegree();
	// preprocess(ego) was called before the start,
	// so pd is for the current network
	const Network * pStart = this->pData()->pNetwork(this->period());
	this->primaryProperties(ego, pStart);
// Now the primary setting is as for the initial network.
//	if (ego < 5) Rprintf("%d %d %d %d  \n", ego,pd , this->lprimDegree);

	int outsidePS = 0; // will be number of ties outside initial primary setting
	for (IncidentTieIterator iter = pNetwork->outTies(ego); // in current network
		iter.valid(); iter.next())
	{
		if (!inPrimarySetting(iter.actor())) // for the initial network!
		{
			outsidePS++;
		}
	}

	if (this->linside)
	{
		if (this->luseSize)
		{
			statistic = outsidePS - pNetwork->outDegree(ego);
		}
		else
		{
			statistic = pNetwork->outDegree(ego) - outsidePS;
		}
	}
	else
	{
		statistic = - outsidePS;
	}

	return statistic;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double PrimaryCompressionEffect::tieStatistic(int alter)
{
	double statistic = 0;
	if (this->linside)
	{
		if (inPrimarySetting(alter))
		{
			if (this->luseSize)
			{
				statistic = - this->llogPrimary;
			}
			else
			{
				statistic = 1;
			}
		}
	}
	else
	{
		if (!inPrimarySetting(alter))
		{
			if (this->luseSize)
			{
				statistic = - this->llogNonPrimary;
			}
			else
			{
				statistic = 1;
			}
		}
	}
	return statistic;
}

}
