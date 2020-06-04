/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateContrastEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateContrastEffect class.
 *****************************************************************************/

#include "CovariateContrastEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */

CovariateContrastEffect::CovariateContrastEffect(
	const EffectInfo * pEffectInfo, const bool plus, const bool minus):
		CovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->lplus = plus;
	this->lminus = minus;
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double CovariateContrastEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;

	double xi = this->covariateValue(actor);
	int deg = this->pNetwork()->outDegree(actor);

	if (this->lplus)
	{
		if (deg > xi)
		{
			contribution = difference * (deg - xi);
		}
	}
	if (this->lminus)
	{
		if (deg < xi)
		{
			contribution = difference * (xi - deg);
		}
	}

	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double CovariateContrastEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;

	if (!this->missingCovariateEitherEnd(ego, this->period()))
	{
		double xi = this->covariateValue(ego);
		int deg = this->pNetwork()->outDegree(ego);

		if (this->lplus)
		{
			if (deg > xi)
			{
				statistic = currentValues[ego] * (deg - xi);
			}
		}
		if (this->lminus)
		{
			if (deg < xi)
			{
				statistic = currentValues[ego] * (xi - deg);
			}
		}
	}

	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double CovariateContrastEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	if ((difference[ego] > 0) && (!this->missingCovariateEitherEnd(ego, this->period())))
	{
		double contrast = this->pNetwork()->outDegree(ego) - this->covariateValue(ego);

		if (this->lplus)
		{
			if (contrast > 0)
			{
				statistic = - difference[ego] * contrast;
			}
		}
		if (this->lminus)
		{
			if (contrast < 0)
			{
				statistic = difference[ego] * contrast;
			}
		}
	}

	return statistic;
}

}
