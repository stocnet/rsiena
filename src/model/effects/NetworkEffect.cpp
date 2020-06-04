/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <R_ext/Print.h>

#include "NetworkEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "model/State.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
NetworkEffect::NetworkEffect(const EffectInfo * pEffectInfo) :
	Effect(pEffectInfo)
{
	this->lpNetwork = 0;
	this->lpNetworkData = 0;
	this->lpNetworkCache = 0;
	this->lstepTypeVal = -1;
	this->lpTwoPathTable = 0;
	this->lpReverseTwoPathTable = 0;
	this->lpInStarTable = 0;
	this->lpOutStarTable = 0;
	this->lpCriticalInStarTable = 0;
	this->lpRRTable = 0;
	this->lpRFTable = 0;
	this->lpRBTable = 0;
	this->lpFRTable = 0;
	this->lpBRTable = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void NetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	Effect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->variableName();

	this->lpNetworkData = pData->pNetworkData(name);

	if (!this->lpNetworkData)
	{
		throw logic_error("Data for network variable '" + name +"' expected.");
	}

	this->lpNetwork = pState->pNetwork(name);
	this->lpNetworkCache = pCache->pNetworkCache(this->lpNetwork);
	this->lstepTypeVal = this->lpNetworkCache->stepTypeValue();

	this->lpTwoPathTable = this->lpNetworkCache->pTwoPathTable();
	this->lpReverseTwoPathTable = this->lpNetworkCache->pReverseTwoPathTable();
	this->lpInStarTable = this->lpNetworkCache->pInStarTable();
	this->lpOutStarTable = this->lpNetworkCache->pOutStarTable();
	this->lpCriticalInStarTable = this->lpNetworkCache->pCriticalInStarTable();
	this->lpRRTable = this->lpNetworkCache->pRRTable();
	this->lpRFTable = this->lpNetworkCache->pRFTable();
	this->lpRBTable = this->lpNetworkCache->pRBTable();
	this->lpFRTable = this->lpNetworkCache->pFRTable();
	this->lpBRTable = this->lpNetworkCache->pBRTable();
}

void NetworkEffect::initialize(const Data * pData, State * pState,
		State * pSimulatedState, int period, Cache * pCache) {
	initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void NetworkEffect::preprocessEgo(int ego)
{
	this->lego = ego;
}


/**
 * Returns if there is a tie from the current ego to the given alter.
 */
bool NetworkEffect::outTieExists(int alter) const
{
	return this->lpNetworkCache->outTieExists(alter);
}


/**
 * Returns if there is a tie from the given alter to the current ego.
 */
bool NetworkEffect::inTieExists(int alter) const
{
	return this->lpNetworkCache->inTieExists(alter);
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function.
 */
double NetworkEffect::evaluationStatistic()
{
	return this->statistic(this->pNetwork());
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function.
 */
pair<double, double * > NetworkEffect::evaluationStatistic(bool needActorStatistics)
{
		return this->statistic(this->pNetwork(), needActorStatistics);
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double NetworkEffect::endowmentStatistic(Network * pLostTieNetwork)
{
	return this->statistic(pLostTieNetwork);
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
pair<double, double *> NetworkEffect::endowmentStatistic(Network * pLostTieNetwork, bool needActorStatistics)
{
	return this->statistic(pLostTieNetwork, needActorStatistics);
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the creation function.
 */
double NetworkEffect::creationStatistic(Network * pGainedTieNetwork)
{
	// We use a trick here. The creation statistic is very similar to the
	// endowment statistic, except that the summation is over the gained
	// ties rather than the lost ties. Normally, you should not override
	// this method in derived classes.

	return this->endowmentStatistic(pGainedTieNetwork);
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the creation function.
 */
pair<double, double *> NetworkEffect::creationStatistic(Network * pGainedTieNetwork, bool needActorStatistics)
{
	return this->endowmentStatistic(pGainedTieNetwork, needActorStatistics);
}



/**
 * Returns if this effect is an ego effect.
 */
bool NetworkEffect::egoEffect() const
{
	// Not an ego effect by default.
	return false;
}


/**
 * A convenience method for implementing statistics for both evaluation and
 * endowment function. It assumes that the statistic can be calculated by
 * iterating over ties (i,j) of a network Y and summing up some terms
 * s_{ij}(X) with respect to another network X, namely,
 * s(X,Y) = sum_{(i,j) \in Y} s_{ij}(X).
 * For evaluation function, X = Y.
 * For endowment function, X is the initial network of the period, and Y is the
 * network of ties that have been lost during the network evolution.
 */
double NetworkEffect::statistic(const Network * pSummationTieNetwork)
{
	return statistic(pSummationTieNetwork, false).first;
}

pair<double,double *> NetworkEffect::statistic(const Network * pSummationTieNetwork, bool needActorStatistics)
{
	this->initializeStatisticCalculation();

	int n = pSummationTieNetwork->n();
	Cache * pCache = this->pCache();
	double statistic = 0;
	double * actorStatistics = 0;
	if(needActorStatistics)
	{
	  actorStatistics = new double[n];
	}
	for (int i = 0; i < n; i++)
	{
		pCache->initialize(i);
		this->preprocessEgo(i);
		this->onNextEgo(i);
		if(needActorStatistics)
		{
			actorStatistics[i] = this->egoStatistic(i, pSummationTieNetwork);
			statistic += actorStatistics[i];
		}
		else
		{
			statistic += this->egoStatistic(i, pSummationTieNetwork);
		}
	}

	this->cleanupStatisticCalculation();
	//Rprintf(" %f sum \n ", statistic);

	return make_pair(statistic,actorStatistics);
}

/**
 * Calculates the statistic corresponding to the given ego. The variable
 * pSummationTieNetwork is the current network in the case of an evaluation
 * effect and the network of lost ties in the case of an endowment effect.
 */
double NetworkEffect::egoStatistic(int ego,
	const Network * pSummationTieNetwork)
{
	double statistic = 0;

	for (IncidentTieIterator iter = pSummationTieNetwork->outTies(ego);
		iter.valid();
		iter.next())
	{
		statistic += this->tieStatistic(iter.actor());
	}
	return statistic;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double NetworkEffect::tieStatistic(int alter)
{
	throw runtime_error("tieStatistic not implemented for " +
		this->pEffectInfo()->effectName());
}


/**
 * This method is called at the start of the calculation of the statistic.
 */
void NetworkEffect::initializeStatisticCalculation()
{
}


/**
 * This method is called at the end of the calculation of the statistic.
 */
void NetworkEffect::cleanupStatisticCalculation()
{
}


/**
 * This method is called right before summing up the contributions of the
 * outgoing ties of the given ego in the calculation of the statistic.
 */
void NetworkEffect::onNextEgo(int ego)
{
}

}
