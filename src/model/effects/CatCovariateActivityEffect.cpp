/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CatCovariateActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CatCovariateActivityEffect class.
 *****************************************************************************/

#include <cmath>
#include "CatCovariateActivityEffect.h"
#include "utils/Utils.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
CatCovariateActivityEffect::CatCovariateActivityEffect(
		const EffectInfo * pEffectInfo) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lpNumberTieValues = 0;
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = (pEffectInfo->internalEffectParameter() >= 2);
}


/**
 * Deallocates this double covariate object.
 */
CatCovariateActivityEffect::~CatCovariateActivityEffect()
{
	delete[] this->lpNumberTieValues;
	this->lpNumberTieValues = 0;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CatCovariateActivityEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);
	int min = int(round(this->covariateMinimum()));
	int max = int(round(this->covariateMaximum()));
	if (min < 0)
	{
		throw logic_error("homXOutAct2: minimum of first covariate is negative");
	}
	if (max > 20)
	{
		throw logic_error("homXOutAct2: first covariate has a maximum which is too large");
	}
	if ((this->pBehaviorData()) || (this->pContinuousData()))
	{
		throw logic_error("homXOutAct2: not implemented for behavioral variables");
	}
	this->lcovMax = max;
	this->lcovMax++;
	this->lpNumberTieValues = new int[this->lcovMax] {};
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling CatCovariateActivityEffect::calculateContribution(...).
 */
void CatCovariateActivityEffect::preprocessEgo(int ego)
{
	CovariateDependentNetworkEffect::preprocessEgo(ego);
	const Network * pNetwork = this->pNetwork();
	for (int b = 0; b < this->lcovMax; b++)
	{
		this->lpNumberTieValues[b]= 0;
	}
	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();
		this->lpNumberTieValues[int(round(this->value(j)))]++;
	}
}

double CatCovariateActivityEffect::changeStat(double d, bool diffSqrt) const
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
 * It is assumed that preprocessEgo(ego) has been called before.
 */
double CatCovariateActivityEffect::calculateContribution(int alter) const
{
	double contribution = 0;
	int altervalue = int(round(this->value(alter)));
	if (altervalue >= 1)
	{
		contribution = this->lpNumberTieValues[altervalue];
		if (this->outTieExists(alter))
		{
			contribution -=2;
		}
		contribution = changeStat(contribution, this->lroot);
	}
	return contribution;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CatCovariateActivityEffect::tieStatistic(int alter)
{
	double contribution = 0;
	if (!(this->missing(alter)))
	{
		int altervalue = int(round(this->value(alter)));
		if (altervalue >= 1)
		{
			contribution = this->lpNumberTieValues[int(round(this->value(alter)))];
		}
		if (this->lroot)
		{
			contribution = this->lsqrtTable->sqrt(contribution);
		}
	}
	return contribution;
}

}
