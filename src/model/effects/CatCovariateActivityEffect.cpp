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
		CatCovariateDependentNetworkEffect(pEffectInfo)
{
	this->lpNumberTieValues = 0;
	this->lsqrtTable = SqrtTable::instance();
	this->lroot = ((int(round(pEffectInfo->internalEffectParameter())) == 2)||
				(int(round(pEffectInfo->internalEffectParameter() == 4))));
	this->laverage = (int(round(pEffectInfo->internalEffectParameter())) >= 3);
}


/**
 * Deallocates this double covariate object.
 */
CatCovariateActivityEffect::~CatCovariateActivityEffect()
{
	delete[] this->lpAllCovariateTies;
	this->lpAllCovariateTies = 0;
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
	CatCovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);
	int min = int(round(this->covariateMinimum()));
	int max = int(round(this->covariateMaximum()));
	if (min < 0)
	{
		throw logic_error("homXOutAct2: minimum of covariate is negative");
	}
	if (max > 20)
	{
		throw logic_error("homXOutAct2: covariate has a maximum which is too large");
	}
	if ((this->pBehaviorData()) || (this->pContinuousData()))
	{
		throw logic_error("homXOutAct2: not implemented for behavioral variables");
	}
	this->lcovMax = max;
	this->lcovMax++;
	this->lpNumberTieValues = new int[this->lcovMax] {};
	this->lpAllCovariateTies  = new double[this->lcovMax] {};
	for (int b = 0; b < this->lcovMax; b++)
	{
		if (this->lroot)
		{
			this->lpAllCovariateTies[b] =
					this->lsqrtTable->sqrt(this->numberCovariateTies(b));
		}
		else
		{
			this->lpAllCovariateTies[b] = this->numberCovariateTies(b);
		}
	}
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling CatCovariateActivityEffect::calculateContribution(...).
 */
void CatCovariateActivityEffect::preprocessEgo(int ego)
{
	CatCovariateDependentNetworkEffect::preprocessEgo(ego);
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
		this->lpNumberTieValues[this->covariateIntValue(j)]++;
	}
}

double CatCovariateActivityEffect::changeStat(double d, bool diffSqrt) const
{
	if (diffSqrt)
	{
		return(((d+1)*this->lsqrtTable->sqrt(d+1)) -
							(d * this->lsqrtTable->sqrt(d)));
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
	int altervalue = this->covariateIntValue(alter);
	if (altervalue >= 1)
	{
		contribution = this->lpNumberTieValues[altervalue];
		if (this->outTieExists(alter))
		{
			contribution -=2;
		}
		contribution = changeStat(contribution, this->lroot);
		if (this->laverage)
		{
			contribution /= this->lpAllCovariateTies[altervalue];
		}
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
		int altervalue = this->covariateIntValue(alter);
		if (altervalue >= 1)
		{
			contribution = this->lpNumberTieValues[altervalue];
			if (this->lroot)
			{
				contribution = this->lsqrtTable->sqrt(contribution);
			}
			if (this->laverage)
			{
				contribution /= this->lpAllCovariateTies[altervalue];
			}
		}
	}
	return contribution;
}

}
