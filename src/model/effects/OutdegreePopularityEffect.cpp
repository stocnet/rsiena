/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreePopularityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutdegreePopularityEffect class.
 *****************************************************************************/

#include <cmath>
#include "data/NetworkLongitudinalData.h"
#include "OutdegreePopularityEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "data/NetworkLongitudinalData.h"
#include "model/EffectInfo.h"
#include "data/Data.h"
#include "utils/Utils.h"


namespace siena
{

/**
 * Constructor.
 * The use of ltrunc is only for clarity in programming.
 */
OutdegreePopularityEffect::OutdegreePopularityEffect(
		const EffectInfo * pEffectInfo, bool root, bool centered, 
					bool threshold, bool trunc) : NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
	this->lcentered = centered;
	this->lcentering = 0.0;
	this->lthreshold = threshold;
	this->ltrunc = trunc;
	if (this->ltrunc)
	{
		int p = pEffectInfo->internalEffectParameter();
		if (this->lroot)
		{
			this->lp = this->lsqrtTable->sqrt(p);
		}
		else
		{
			this->lp = p;
		}
	}
	else
	{
		this->lp = 0;
	}
	this->lvariableName = pEffectInfo->variableName();
// centering and root cannot occur simultaneously
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void OutdegreePopularityEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	if (this->lcentered)
	{
		NetworkLongitudinalData * pNetworkData =
				pData->pNetworkData(this->lvariableName);
		this->lcentering = pNetworkData->averageOutDegree();
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutdegreePopularityEffect::calculateContribution(int alter) const
{
	int degree = this->pNetwork()->outDegree(alter);
	double change = degree;

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(degree);
	}
		
	if (this->ltrunc)
	{
		if (change <= this->lp + EPSILON)
		{
			change = 0;
		}
		else
		{
			change = change - this->lp;
		}
		if ((this->lthreshold) && (change > 0.0001))
		{
			change = 1;
		}
	}
	else	
	{
		change = change - this->lcentering;
	}
	
	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double OutdegreePopularityEffect::tieStatistic(int alter)
{
	return this->calculateContribution(alter);
}

}
