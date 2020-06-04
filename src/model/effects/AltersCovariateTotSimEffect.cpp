/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateTotSimEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AltersCovariateTotSimEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <cstdlib>

#include "AltersCovariateTotSimEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

namespace siena {

/**
 * Constructor.
 */
AltersCovariateTotSimEffect::AltersCovariateTotSimEffect(
		const EffectInfo * pEffectInfo) :
	CovariateAndNetworkBehaviorEffect(pEffectInfo) {
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AltersCovariateTotSimEffect::calculateChangeContribution(int actor,
		int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0) // otherwise, nothing to calculate...
	{

		int oldValue = this->value(actor); // ego's behavior value before moving on behavior scale
		int newValue = oldValue + difference; // ego's behavior value after moving on behavior scale
		double totalChange = 0;                   // will keep track of changes

		for (IncidentTieIterator iter = pNetwork->outTies(actor); // loops over outgoing ties of ego
		iter.valid(); iter.next()) {
			int j = iter.actor();                // identifies alter
			int alterValue = this->value(j); // identifies behavior value of alter
			int change = std::abs(oldValue - alterValue)
					- std::abs(newValue - alterValue);
			// calculates impact of ego's movement on absolute difference to alter

			totalChange += change * this->covariateValue(j); // weighting change statistic acc. to covariate value of alter
		}

		contribution = ((double) totalChange) / this->range(); // standardize by behavior range for SIMILARITY function

	}

	return contribution;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AltersCovariateTotSimEffect::egoStatistic(int ego,
		double * currentValues)
{
	const Network * pNetwork = this->pNetwork();
	double statistic = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
			iter.valid();
			iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j)
				&& !this->missing(this->period() + 1, j)
				&& !this->missingCovariate(j, this->period())) {
			double tieStatistic = this->similarity(currentValues[ego],
					currentValues[j]);

			statistic += tieStatistic * this->covariateValue(j);
		}
	}

	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AltersCovariateTotSimEffect::egoEndowmentStatistic(int ego,
		const int * difference, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	if (difference[ego] > 0 && !this->missingDummy(ego)
			&& (pNetwork->outDegree(ego) > 0)) // otherwise, nothing to calculate...
	{

		int oldValue = this->value(ego); // ego's behavior value before moving on behavior scale
		int newValue = oldValue + difference[ego]; // ego's behavior value after moving on behavior scale
		double totalChange = 0;                   // will keep track of changes

		for (IncidentTieIterator iter = pNetwork->outTies(ego); // loops over outgoing ties of ego
				iter.valid(); iter.next()) {
			int j = iter.actor();                // identifies alter
			int alterValue = this->value(j); // identifies behavior value of alter
			int change = std::abs(oldValue - alterValue)
				- std::abs(newValue - alterValue);
			// calculates impact of ego's movement on absolute difference to alter

			totalChange += change * this->covariateValue(j); // weighting change statistic acc. to covariate value of alter
		}

		statistic -= totalChange / this->range(); // standardize by behavior range for SIMILARITY function

	}

	return statistic;
}

}
