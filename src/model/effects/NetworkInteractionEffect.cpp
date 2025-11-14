/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkInteractionEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * NetworkInteractionEffect.
 *****************************************************************************/

#include "NetworkInteractionEffect.h"

namespace siena
{

/**
 * Constructs a new interaction effect between the given effects.
 * The parameter pEffect3 should be 0 for two-way interactions.
 * This effect takes the ownership of the given effects, which mean
 * that the given effects will be destroyed as soon as this
 * effect is destroyed.
 */
NetworkInteractionEffect::NetworkInteractionEffect(
	const EffectInfo * pEffectInfo,
	NetworkEffect * pEffect1,
	NetworkEffect * pEffect2,
	NetworkEffect * pEffect3) : NetworkEffect(pEffectInfo)
{
	this->lpEffect1 = pEffect1;
	this->lpEffect2 = pEffect2;
	this->lpEffect3 = pEffect3;
}


/**
 * Deallocates this effects.
 */
NetworkInteractionEffect::~NetworkInteractionEffect()
{
	delete this->lpEffect1;
	delete this->lpEffect2;
	delete this->lpEffect3;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void NetworkInteractionEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	this->lpEffect1->initialize(pData, pState, period, pCache);
	this->lpEffect2->initialize(pData, pState, period, pCache);

	if (this->lpEffect3)
	{
		this->lpEffect3->initialize(pData, pState, period, pCache);
	}
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void NetworkInteractionEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);

	this->lpEffect1->preprocessEgo(ego);
	this->lpEffect2->preprocessEgo(ego);

	if (this->lpEffect3)
	{
		this->lpEffect3->preprocessEgo(ego);
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double NetworkInteractionEffect::calculateContribution(int alter) const
{
	double contribution =
		this->lpEffect1->calculateContribution(alter) *
		this->lpEffect2->calculateContribution(alter);

	if (this->lpEffect3)
	{
		contribution *= this->lpEffect3->calculateContribution(alter);
	}

	return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego. The variable
 * pSummationTieNetwork is the current network in the case of an evaluation
 * effect and the network of lost ties in the case of an endowment effect.
 */
double NetworkInteractionEffect::egoStatistic(int ego,
	const Network * pSummationTieNetwork)
{
	int egoEffectCount = 0;
	int effectCount = 2;

	if (this->lpEffect1->egoEffect())
	{
		egoEffectCount++;
	}

	if (this->lpEffect2->egoEffect())
	{
		egoEffectCount++;
	}

	if (this->lpEffect3)
	{
		effectCount++;

		if (this->lpEffect3->egoEffect())
		{
			egoEffectCount++;
		}
	}

	double statistic;

	if (egoEffectCount == effectCount - 1)
	{
		// Special case. The only effect that is not an ego effect can be
		// an arbitrary network effect.
		// Distinguishing this special case leads only to a 
		// minor gain in efficiency.

		// In what follows, we use the fact that tieStatistic(alter) is
		// the same for an ego effect regardless of the alter.

		if (this->lpEffect1->egoEffect())
		{
			statistic = this->lpEffect1->tieStatistic(0);
		}
		else
		{
			statistic =
				this->lpEffect1->egoStatistic(ego, pSummationTieNetwork);
		}

		if (this->lpEffect2->egoEffect())
		{
			statistic *= this->lpEffect2->tieStatistic(0);
		}
		else
		{
			statistic *=
				this->lpEffect2->egoStatistic(ego, pSummationTieNetwork);
		}

		if (this->lpEffect3)
		{
			if (this->lpEffect3->egoEffect())
			{
				statistic *= this->lpEffect3->tieStatistic(0);
			}
			else
			{
				statistic *=
					this->lpEffect3->egoStatistic(ego, pSummationTieNetwork);
			}
		}
	}
	else
	{
		// The regular case: All interacting effects are assumed to be dyadic
		statistic = NetworkEffect::egoStatistic(ego, pSummationTieNetwork);
		// this means that the sum of
		// NetworkInteractionEffect::tieStatistic (below!) will be used
	}

	return statistic;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double NetworkInteractionEffect::tieStatistic(int alter)
{
	double statistic =
		this->lpEffect1->tieStatistic(alter) *
		this->lpEffect2->tieStatistic(alter);

	if (this->lpEffect3)
	{
		statistic *= this->lpEffect3->tieStatistic(alter);
	}

	return statistic;
}


/**
 * This method is called at the start of the calculation of the statistic.
 */
void NetworkInteractionEffect::initializeStatisticCalculation()
{
	this->lpEffect1->initializeStatisticCalculation();
	this->lpEffect2->initializeStatisticCalculation();

	if (this->lpEffect3)
	{
		this->lpEffect3->initializeStatisticCalculation();
	}
}


/**
 * This method is called at the end of the calculation of the statistic.
 */
void NetworkInteractionEffect::cleanupStatisticCalculation()
{
	this->lpEffect1->cleanupStatisticCalculation();
	this->lpEffect2->cleanupStatisticCalculation();

	if (this->lpEffect3)
	{
		this->lpEffect3->cleanupStatisticCalculation();
	}
}


/**
 * This method is called right before summing up the contributions of the
 * outgoing ties of the given ego in the calculation of the statistic.
 */
void NetworkInteractionEffect::onNextEgo(int ego)
{
	this->lpEffect1->onNextEgo(ego);
	this->lpEffect2->onNextEgo(ego);

	if (this->lpEffect3)
	{
		this->lpEffect3->onNextEgo(ego);
	}
}


/**
 * Returns if this effect is an ego effect.
 */
bool NetworkInteractionEffect::egoEffect() const
{
	bool rc = this->lpEffect1->egoEffect() && this->lpEffect2->egoEffect();

	if (this->lpEffect3)
	{
		rc = rc && this->lpEffect3->egoEffect();
	}

	return true;
}

}
