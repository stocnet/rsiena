/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAvSimEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AltersCovariateAvAltEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>

#include "AltersCovariateAvAltEffect.h"
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
AltersCovariateAvAltEffect::AltersCovariateAvAltEffect(
	const EffectInfo * pEffectInfo, bool divide, bool same) :
	CovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
	this->lsame = same;
	// Indicates whether there will be restriction to same covariate
}

/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void AltersCovariateAvAltEffect::preprocessEgo(int ego)
{
// do not repeat all this:	CovariateAndNetworkBehaviorEffect::preprocessEgo(ego);
    this->lTotalAlterValue = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(ego) > 0)
	{
		double totalAlterValue = 0;
		double covEgo = this->covariateValue(ego);
		int neighborCount = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(ego);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();                // identifies alter
			if (this->lsame)
			{
				if (fabs(this->covariateValue(j) - covEgo) < 1e-6)
				{
					totalAlterValue +=  this->centeredValue(j);
					neighborCount++;
				}
			}
			else
			{
				totalAlterValue +=  this->centeredValue(j) * this->covariateValue(j);
				neighborCount++;
			}
		}
		if  (neighborCount > 0)
		{
			if (this->ldivide)
			{
				this->lTotalAlterValue = totalAlterValue / neighborCount;
			}
			else
			{
			this->lTotalAlterValue = totalAlterValue;
			}
		}
	}	
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AltersCovariateAvAltEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->lTotalAlterValue;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AltersCovariateAvAltEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;
	double covEgo = this->covariateValue(ego);

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j) &&
			!this->missingCovariate(j,this->period()))
		{	
			if (this->lsame)
			{
				if (fabs(this->covariateValue(j) - covEgo) < 1e-6)
				{
					statistic += currentValues[j];
					neighborCount++;
				}
			}
			else
			{
				statistic += currentValues[j] * this->covariateValue(j);
				neighborCount++;
			}			
		}
	}

	if (neighborCount > 0) 
	{
		if (this->ldivide)
		{
			statistic *= currentValues[ego] / neighborCount;
		}
		else
		{
			statistic *= currentValues[ego];
		}
	}
	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AltersCovariateAvAltEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;
	int neighborCount = 0;
	double covEgo = this->covariateValue(ego);
	const Network * pNetwork = this->pNetwork();

	if ((difference[ego] > 0) && !this->missingDummy(ego) && (pNetwork->outDegree(ego) > 0)) // otherwise, nothing to calculate...
	{
		double totalAlterValue = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(ego);
				iter.valid();
				iter.next())
		{
			int j = iter.actor();                // identifies alter
			if (this->lsame)
			{
				if (fabs(this->covariateValue(j) - covEgo) < 1e-6)
				{
					totalAlterValue += currentValues[j];
					neighborCount++;
				}
			}
			else
			{
				totalAlterValue += currentValues[j] * this->covariateValue(j);
				neighborCount++;
			}			
		}

		if (neighborCount > 0)
		{
			if (this->ldivide)
			{
				statistic = - difference[ego] * totalAlterValue / neighborCount;
			}
			else
			{
				statistic = - difference[ego] * totalAlterValue;
			}
		}
	}
	return statistic;
}


}
