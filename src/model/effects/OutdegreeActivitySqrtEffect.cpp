/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreeActivitySqrtEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutdegreeActivitySqrtEffect class.
 *****************************************************************************/

#include <cmath>

#include "OutdegreeActivitySqrtEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "data/NetworkLongitudinalData.h"

namespace siena
{

/**
 * Constructor.
 */
OutdegreeActivitySqrtEffect::OutdegreeActivitySqrtEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutdegreeActivitySqrtEffect::calculateContribution(int alter)
	const
{
	// Current out-degree

	int degree = this->pNetwork()->outDegree(this->ego());

	if (this->outTieExists(alter))
	{
		// The tie is going to be withdrawn, so the degree decreases

		return degree * this->lsqrtTable->sqrt(degree) -
			(degree - 1) * this->lsqrtTable->sqrt(degree - 1);
	}
	else
	{
		// The tie is going to be introduced, so the degree increases

		return (degree + 1) * this->lsqrtTable->sqrt(degree + 1) -
			degree * this->lsqrtTable->sqrt(degree);
	}
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double OutdegreeActivitySqrtEffect::tieStatistic(int alter)
{
	return this->lsqrtTable->sqrt(this->pNetwork()->outDegree(this->ego()));
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double OutdegreeActivitySqrtEffect::endowmentStatistic(Network * pLostTieNetwork)
{
	double statistic = 0;
	const Network* pStart = this->pData()->pNetwork(this->period());
	int n = pStart->n();
	for (int i = 0; i < n; i++)
	{
		double routdeg =  this->lsqrtTable->sqrt(pStart->outDegree(i));
		statistic += routdeg * pLostTieNetwork->outDegree(i);
	}
	return statistic;
}

}
