/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PrimarySettingEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * PrimarySettingEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include <string>
#include <R_ext/Print.h>
#include "NetworkEffect.h"
#include "PrimarySettingEffect.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "network/UnionNeighborIterator.h"
#include "data/NetworkLongitudinalData.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
PrimarySettingEffect::PrimarySettingEffect(const EffectInfo * pEffectInfo,
		bool difference, bool logar, bool root, bool inv): NetworkWithPrimaryEffect(pEffectInfo)
{
	this->lparameter = pEffectInfo->internalEffectParameter();
	this->ldifference = difference;
	this->llogar = logar;
	this->lroot = root;
	this->linv = inv;
	this->levalBoth = (this->lparameter >= 2);
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void PrimarySettingEffect::initialize(const Data * pData,
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
void PrimarySettingEffect::preprocessEgo(int ego)
{
	NetworkWithPrimaryEffect::preprocessEgo(ego);
}

/**
 * Returns if this effect is an ego effect.
 */
bool PrimarySettingEffect::egoEffect() const
{
	return true;
}


double PrimarySettingEffect::transform(int value) const
/**
 * Transforms value
 */
{
	if (value < 0)
	{
		throw logic_error("negative value in PrimarySettingEffect::transform: value= " +
			to_string(value));
	}
	double contribution = value;
	if (this->llogar)
	{
		contribution = log(contribution+1);
	}
	else if (this->lroot)
	{
		contribution = sqrt(contribution);
	}
	else if (this->linv)
	{
		contribution = 1/(contribution + 1);
	}
	return contribution;
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double PrimarySettingEffect::calculateContribution(int alter) const
{
	double contribution = 0;
	int contr = this->primaryDegree();
	int outdeg = this->pNetwork()->outDegree(this->ego());
	if (this->ldifference)
	{
		contr -= outdeg;
	}
	contribution = this->transform(contr);
	return contribution;
}

/**
 * The contribution of ego to the statistic.
 * outdegree multiplied by primary setting degree;
 * if evalBoth, the sum of the primary setting degree in the current and the starting network;
 * if !evalBoth, only in the starting network.
 */
double PrimarySettingEffect::egoStatistic(int ego, const Network * pNetwork)
{
	double statistic = 0;
	int pd = this->primaryDegree();
	// preprocess(ego) was called before the start,
	// so pd is for the current network;
	// taken now because this->primaryDegree will be changed
	const Network * pStart = this->pData()->pNetwork(this->period());
	this->primaryProperties(ego, pStart);
//	if (ego < 5) Rprintf("%d %d %d %d  \n", ego,pd , this->primaryDegree());

	int val = this->primaryDegree(); // this one is for initial network
	int addIfBoth = 0;
	double raddIfBoth = 0.0;
	if (this->ldifference)
	{
		val -= pStart->outDegree(ego); // minus outdegree for initial network
	}
	if (this->levalBoth)
	{
		addIfBoth = pd;
		if (this->ldifference)
		{
			addIfBoth -= this->pNetwork()->outDegree(ego); // minus outdegree for current network
		}
		raddIfBoth = this->transform(addIfBoth);
	}
	statistic = (this->transform(val) + raddIfBoth) * pNetwork->outDegree(ego);
// note that this->pNetwork() is the current network, whereas
// pNetwork is the network in the call of egoStatistic,
// which will be the difference network in the case of endowment or creation.
	return statistic;
}

/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double PrimarySettingEffect::tieStatistic(int alter)
{
	// evalBoth not implemented here.
	// It should really be done through a component of OneModeNetworkLongitudinalData.
	int pd = this->primaryDegree();
	if (this->ldifference)
	{
		pd -= this->pNetwork()->outDegree(this->ego());
	}
	return this->transform(pd);
}

}
