/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateAvAltEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DyadicCovariateAvAltEffect class.
 *****************************************************************************/

#include <cmath>
#include "DyadicCovariateAvAltEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
DyadicCovariateAvAltEffect::DyadicCovariateAvAltEffect(
	const EffectInfo * pEffectInfo, bool divide, bool asWeight, bool outgoing) :
	DyadicCovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the degree
	this->lasWeight = asWeight;
	// Indicates that the dyadic covariate is used as a weight;
	// if not, used as the variable.
	this->lpar2 = (pEffectInfo->internalEffectParameter() >= 2);
	// Indicates wether in or outgoing ties should be used
	this->loutgoing = outgoing;
	// specifies type of denominator
	if (!lasWeight) {lpar2 = false;}
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double DyadicCovariateAvAltEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();
	IncidentTieIterator iter;
	bool posDegree = true;
	if (loutgoing)
	{
		iter = pNetwork->outTies(actor);
		posDegree = (pNetwork->outDegree(actor) > 0);
	}
	else
	{
		iter = pNetwork->inTies(actor);		
		posDegree = (pNetwork->inDegree(actor) > 0);
	}
	
	if (posDegree)
	{
		double totalAlterValue = 0;
		double totalWeightValue = 0;
	
		if (this->ldivide) // av
		{
			for ( ; iter.valid(); iter.next())
			{
					int j = iter.actor();                // identifies alter
					double dycova = this->dycoValue(actor, j);
					if (lasWeight)
					{
						totalAlterValue += (double) this->centeredValue(j) * dycova;
					}
					else
					{
						totalAlterValue += (double) dycova;
					}
					if (lpar2)
					{
						totalWeightValue += (double) dycova;
					}
					else
					{
						totalWeightValue += 1;
					}
				}
				if (fabs(totalWeightValue) > EPSILON)  //  normally this will be a comparison of 0 against >= 1
				{
					contribution = (double) difference * totalAlterValue / totalWeightValue;
				}
			}
			else // tot
			{
			for ( ;	iter.valid(); iter.next())
				{
					int j = iter.actor();                // identifies alter
					double dycova = this->dycoValue(actor, j);
					if (lasWeight)
					{
						totalAlterValue += (this->centeredValue(j) * dycova);
					}
					else
					{
						totalAlterValue += (double) dycova;
					}
				}
		contribution = (double) difference * totalAlterValue;
		}
	}
	return contribution;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double DyadicCovariateAvAltEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	double totalWeightValue = 0;
	IncidentTieIterator iter;
	if (loutgoing)
	{
		iter = pNetwork->outTies(ego);
	}
	else
	{
		iter = pNetwork->inTies(ego);		
	}
	
	for ( ; iter.valid(); iter.next())
		{
			int j = iter.actor();
			if (!this->missingDyCo(ego,j))
			{
				double dycova = this->dycoValue(ego, j);
				if (lasWeight)
				{
					statistic += currentValues[j] * dycova;
				}
				else
				{
					statistic += (double) dycova;
				}
				if (lpar2)
				{
					totalWeightValue += (double) dycova;
				}
				else
				{
					totalWeightValue += 1;
				}
			}
		}
	
	statistic *= currentValues[ego];
	if ((this->ldivide) && (fabs(totalWeightValue) > EPSILON)) // normally, comparison between 0 and >= 1
	{
		statistic /= totalWeightValue;
	}

	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double DyadicCovariateAvAltEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	double totalWeightValue = 0;
	IncidentTieIterator iter;
	if (loutgoing)
	{
		iter = pNetwork->outTies(ego);
	}
	else
	{
		iter = pNetwork->inTies(ego);		
	}
	
	if (difference[ego] > 0 && (pNetwork->outDegree(ego) > 0)) // otherwise, nothing to calculate...
	{
		double thisStatistic = 0;
		double previousStatistic = 0;
		for ( ;	iter.valid(); iter.next())
			{
				int j = iter.actor();
				if (!this->missingDyCo(ego,j))
				{
					double dycova = this->dycoValue(ego, j);
					double alterValue = dycova * currentValues[j];
					double alterPreviousValue = dycova *
									(currentValues[j] + difference[j]);
					thisStatistic += alterValue;
					previousStatistic += alterPreviousValue;
					if (lpar2)
					{
						totalWeightValue += (double) dycova;
					}
					else
					{
						totalWeightValue += 1;
					}
				}
			}
			thisStatistic *= currentValues[ego];
			previousStatistic *= (currentValues[ego] + difference[ego]);
			statistic = (thisStatistic - previousStatistic) * currentValues[ego];
			if ((this->ldivide) && (fabs(totalWeightValue) > 1e-15))// normally, comparison between 0 and >= 1
			{
				statistic /= totalWeightValue;
		}
	}
	
	return statistic;
}

}
