/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageGroupEgoEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageGroupEgoEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <string>

#include "AverageGroupEgoEffect.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/Data.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
AverageGroupEgoEffect::AverageGroupEgoEffect(const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
//	this->lcenterMean = (pEffectInfo->internalEffectParameter() <= 0.5);
//	if (!this->lcenterMean)
//	{
//		this->lcenteringValue = pEffectInfo->internalEffectParameter();
//	}
//	else
//	{
//		this->lcenteringValue = 0.0;
//	}
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void AverageGroupEgoEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->interactionName1();
	this->lpBehaviorData = pData->pBehaviorData(name);
	this->lnm = this->pNetwork()->m();
	this->lGroupMean = 0.0;
	this->lperiod = period;

	if (this->pConstantCovariate())
	{
//		this->lGroupMean = this->pConstantCovariate()->mean();
		throw logic_error(
				"avGroupEgoX effect not meaningful for constant covariate '" + name + "'.");
	}
	else if (this->pChangingCovariate())
	{
		this->lGroupMean = 0.0;
		int mnonmis = 0;
		for (int i = 0; i < this->lnm; i++)
		{
			if (!this->pChangingCovariate()->missing(i, this->lperiod))
			{
				this->lGroupMean += this->pChangingCovariate()->value(i, this->lperiod);
				mnonmis++;
			}
		}
		if (mnonmis > 0)
		{
			this->lGroupMean /= mnonmis;
		}
	}
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void AverageGroupEgoEffect::preprocessEgo(int ego)
{
	CovariateDependentNetworkEffect::preprocessEgo(ego);
	if (this->lpBehaviorData)
	{
		this->lGroupMean = 0.0;
		for (int i = 0; i < this->lnm; i++)
		{
			this->lGroupMean += this->value(i);
// according to CovariateDependentNetworkEffect.cpp this is a value
// centered by this->lpBehaviorData->overallMean()
		}
		this->lGroupMean /= this->lnm;
	}
}



/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double AverageGroupEgoEffect::calculateContribution(int alter) const
{
	return this->lGroupMean;
}


/**
 * Returns if this effect is an ego effect.
 */
bool AverageGroupEgoEffect::egoEffect() const
{
	return true;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double AverageGroupEgoEffect::tieStatistic(int alter)
{
	return this->lGroupMean;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double AverageGroupEgoEffect::egoStatistic(int ego,	const Network * pNetwork)
{
	return pNetwork->outDegree(ego) * this->lGroupMean;
}

}

