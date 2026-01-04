/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAverageEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AltersCovariateAverageEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>

#include "AltersCovariateAverageEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AltersCovariateAverageEffect::AltersCovariateAverageEffect(
	const EffectInfo * pEffectInfo, bool divide, bool same, bool outgoing) :
	CovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// ldivide indicates whether there will be division by the outdegree of ego
	this->lsame = same;
	// same indicates whether the influence will be from same-X alters
	this->loutgoing = outgoing;
	// outgoing indicates whether the influence will be from out-alters
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void AltersCovariateAverageEffect::preprocessEgo(int ego)
{
// do not repeat all this:	CovariateAndNetworkBehaviorEffect::preprocessEgo(ego);
    this->lTotalAlterValue = 0;
	const Network * pNetwork = this->pNetwork();	
	int neighborCount = 0;
	
	IncidentTieIterator iter;
	if (this->loutgoing)
	{
		iter = pNetwork->outTies(ego);
	}
	else
	{
		iter = pNetwork->inTies(ego);		
	}

	if (pNetwork->outDegree(ego) > 0)
	{
		double covEgo = this->covariateValue(ego);
		for ( ; iter.valid(); iter.next())
		{
			int j = iter.actor();                // identifies alter
			double alterValue = 0;
			if (this->lsame)
			{
				if (fabs(this->covariateValue(j) - covEgo) < 1e-6)
				{
					alterValue++;				
				}
			}
			else
			{
				alterValue = this->covariateValue(j);
				
			}
			neighborCount++;
			this->lTotalAlterValue += alterValue;
		}

		if ((this->ldivide) & (neighborCount > 0))
		{
			this->lTotalAlterValue /= neighborCount;
		}
	}
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AltersCovariateAverageEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->lTotalAlterValue;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AltersCovariateAverageEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;
	bool misEgo = (this->missing(this->period() + 1, ego) || 				
						this->missing(this->period() + 1, ego));
	double covEgo = this->covariateValue(ego);

	IncidentTieIterator iter;
	if (this->loutgoing)
	{
		iter = pNetwork->outTies(ego);
	}
	else
	{
		iter = pNetwork->inTies(ego);		
	}
	
	if (!(this->lsame && misEgo))
	{
		for ( ;  iter.valid(); iter.next())
		{
			int j = iter.actor();

			if (!this->missing(this->period(), j) &&
				!this->missing(this->period() + 1, j) &&
				!this->missingCovariate(j,this->period()))
			{
				if (this->lsame)
				{
					if ((fabs(this->covariateValue(j) - covEgo) < 1e-6) && (!misEgo))
					{
						statistic ++;			
					}
				}
				else
				{
					statistic += this->covariateValue(j);
				}
				neighborCount++;	
			}
		}
		statistic *= currentValues[ego];
	}

	if ((neighborCount > 0) && (this->ldivide))
	{
		statistic /= neighborCount;
	}

	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AltersCovariateAverageEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int degreeEgo = 0;
	if (this->loutgoing)
	{
		degreeEgo = pNetwork->outDegree(ego);
	}
	else
	{
		degreeEgo = pNetwork->inDegree(ego);	
	}

	if (difference[ego] > 0 && !this->missingDummy(ego) && (degreeEgo > 0)) // otherwise, nothing to calculate...
	{
		int neighborCount = 0;
		double covEgo = this->covariateValue(ego);

		IncidentTieIterator iter;
		if (this->loutgoing)
		{
			iter = pNetwork->outTies(ego);
		}
		else
		{
			iter = pNetwork->inTies(ego);		
		}
		
		for ( ; iter.valid(); iter.next())
		{
			int j = iter.actor();                // identifies alter
			if (this->lsame)
			{
				if (fabs(this->covariateValue(j) - covEgo) < 1e-6)
				{
					statistic ++;	
				}
			}
			else
			{
				statistic += this->covariateValue(j);
			}		
			neighborCount++;	
		}
		
		statistic *= difference[ego];

		if (this->ldivide)
		{
			statistic /= neighborCount;
		}
	}
	return -statistic;
}


}
