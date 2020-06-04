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
	const EffectInfo * pEffectInfo, bool divide) :
	CovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AltersCovariateAvAltEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{

		double totalAlterValue = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();                // identifies alter
			double alterValue = this->centeredValue(j) * this->covariateValue(j);
			totalAlterValue += alterValue;
		}

		if (this->ldivide)
		{
			contribution = difference * totalAlterValue /
				pNetwork->outDegree(actor);
		}
		else
		{
			contribution = difference * totalAlterValue;
		}
	}

	return contribution;
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

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j) &&
			!this->missingCovariate(j,this->period()))
		{
			statistic += currentValues[j] * this->covariateValue(j);
			neighborCount++;
		}
	}

	if ((neighborCount > 0) && (this->ldivide))
	{
		statistic *= currentValues[ego] / neighborCount;
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
	const Network * pNetwork = this->pNetwork();

	if (difference[ego] > 0 && !this->missingDummy(ego) && (pNetwork->outDegree(ego) > 0)) // otherwise, nothing to calculate...
	{
		double totalAlterValue = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(ego);
				iter.valid();
				iter.next())
		{
			int j = iter.actor();                // identifies alter
			double alterValue = this->centeredValue(j) * this->covariateValue(j);
			totalAlterValue += alterValue;
		}

		if (this->ldivide)
		{
			statistic -= difference[ego] * totalAlterValue /
				pNetwork->outDegree(ego);
		}
		else
		{
			statistic -= difference[ego] * totalAlterValue;
		}
	}
	return statistic;
}


}
