/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CrossCovariateActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CrossCovariateActivityEffect class.
 *****************************************************************************/

#include <cmath>
#include "CrossCovariateActivityEffect.h"
#include "utils/SqrtTable.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"


namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
CrossCovariateActivityEffect::CrossCovariateActivityEffect(
		const EffectInfo * pEffectInfo, bool recip) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lrecip = recip;
	this->lsqrt = ((pEffectInfo->internalEffectParameter() == 2) ||
						(pEffectInfo->internalEffectParameter() == 3));
	this->lsqrtTable = SqrtTable::instance();
	this->lthree = (pEffectInfo->internalEffectParameter() == 3);
}


/**
 * Determines the condition combining everything
 */

bool CrossCovariateActivityEffect::lcondition1(int theAlter, double theOwnValue) const
{
	return ((fabs(this->value(theAlter) - theOwnValue) < EPSILON) &&
			(!lrecip || (this->inTieExists(theAlter))));
}

double CrossCovariateActivityEffect::changeStat(double d, bool diffSqrt) const
{
	if (diffSqrt)
	{
		return(((d+1)*this->lsqrtTable->sqrt(d+1)) - (d * this->lsqrtTable->sqrt(d)));
	}
	else
	{
		return((2*d) + 1);
	}	
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
// For recip this is not very efficient, because each alter has the same contribution;
// could be made more efficient by preprocessEgo?
double CrossCovariateActivityEffect::calculateContribution(int alter) const
{
	double myvalue = this->value(this->ego());
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();		
	int degsame = 0;
	int degdiff = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
		iter.valid();
		iter.next())
	{
		// Get the receiver of the outgoing tie.
		int h = iter.actor();
		if (lcondition1(h, myvalue)) // (myvalue == this->value(h))
		{
			degsame++;
		}
		else
		{
			degdiff++;
		}
	}
	if (this->outTieExists(alter))
	{
		if (lcondition1(alter, myvalue)) // (myvalue == this->value(alter))
		{
			degsame--;
		}
		else
		{
			degdiff--;
		}
	}

	if (this->lthree)
	{
			if (lcondition1(alter, myvalue)) // (myvalue == this->value(alter))
			{
				contribution = changeStat(degsame, true) * degdiff;			
			}
			else
			{
				contribution = changeStat(degdiff, false) *  this->lsqrtTable->sqrt(degsame);		
			}			
	}
	else
	{
		if (this->lsqrt)
		{
			if (lcondition1(alter, myvalue)) // (myvalue == this->value(alter))
			{
				contribution = changeStat(degsame, false) * this->lsqrtTable->sqrt(degdiff);			
			}
			else
			{
				contribution = changeStat(degdiff, true) *  degsame;		
			}			
		}
		else
		{
			if (lcondition1(alter, myvalue)) // (myvalue == this->value(alter))
			{
				contribution = changeStat(degsame, false) * degdiff;			
			}
			else
			{
				contribution = changeStat(degdiff, false) * degsame;		
			}			
		}
	}
	return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego. The parameter
 * pNetwork is always the current network as there are no endowment effects
 * of this kind.
 */
double CrossCovariateActivityEffect::egoStatistic(int ego,
	const Network * pNetwork)
{
	double statistic = 0;	
	int degsame = 0;
	int degdiff = 0;
	double myvalue = this->value(ego);

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		iter.valid();
		iter.next())
	{
		// Get the receiver of the outgoing tie.
		int h = iter.actor();
		if (lcondition1(h, myvalue)) // (myvalue == this->value(h))
		{
			degsame++;
		}
		else
		{
			degdiff++;
		}
	}	
	if (this->lthree)
	{
		statistic = this->lsqrtTable->sqrt(degsame) * degdiff;
	}
	else if (this->lsqrt)
	{
		statistic = this->lsqrtTable->sqrt(degdiff) * degsame;		
	}
	else
	{
		statistic = degsame * degdiff;
	}
	return statistic;
}

/**
 * Returns if this effect is an ego effect.
 */
bool CrossCovariateActivityEffect::egoEffect() const
{
	return true;
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double CrossCovariateActivityEffect::endowmentStatistic(Network * pLostTieNetwork)
{
	throw std::logic_error(
		"outdegree activity differen-same covariate effect: Endowment effect not supported.");
}


}
