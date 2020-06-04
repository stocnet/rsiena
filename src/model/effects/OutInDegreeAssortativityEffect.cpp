/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutInDegreeAssortativityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutInDegreeAssortativityEffect class.
 *****************************************************************************/

#include <cmath>

#include "OutInDegreeAssortativityEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "network/TieIterator.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
OutInDegreeAssortativityEffect::OutInDegreeAssortativityEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lroot = pEffectInfo->internalEffectParameter() == 2;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego.
 */
void OutInDegreeAssortativityEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);

	const Network * pNetwork = this->pNetwork();
	this->ldegree = pNetwork->outDegree(ego);

	if (this->lroot)
	{
		this->lsqrtDegree = this->lsqrtTable->sqrt(this->ldegree);
		this->lsqrtDegreePlus = this->lsqrtTable->sqrt(this->ldegree + 1);

		if (this->ldegree >= 1)
		{
			this->lsqrtDegreeMinus = this->lsqrtTable->sqrt(this->ldegree - 1);
		}
	}

	this->lneighborDegreeSum = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		iter.valid();
		iter.next())
	{
		int alterDegree = pNetwork->inDegree(iter.actor());

		if (this->lroot)
		{
			this->lneighborDegreeSum += this->lsqrtTable->sqrt(alterDegree);
		}
		else
		{
			this->lneighborDegreeSum += alterDegree;
		}
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutInDegreeAssortativityEffect::calculateContribution(int alter)
	const
{
	double change = 0;
	int alterDegree = this->pNetwork()->inDegree(alter);

	if (this->outTieExists(alter))
	{
		if (this->lroot)
		{
			double sqrtAlterDegree = this->lsqrtTable->sqrt(alterDegree);

			change =
				(this->lneighborDegreeSum - sqrtAlterDegree) *
					(this->lsqrtDegree - this->lsqrtDegreeMinus) +
				this->lsqrtDegree * sqrtAlterDegree;
		}
		else
		{
			change =
				(this->lneighborDegreeSum - alterDegree) +
				this->ldegree * alterDegree;
		}
	}
	else
	{
		if (this->lroot)
		{
			change =
				this->lneighborDegreeSum *
					(this->lsqrtDegreePlus - this->lsqrtDegree) +
				this->lsqrtDegreePlus *
					this->lsqrtTable->sqrt(alterDegree + 1);
		}
		else
		{
			change =
				this->lneighborDegreeSum +
				(this->ldegree + 1) * (alterDegree + 1);
		}
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double OutInDegreeAssortativityEffect::tieStatistic(int alter)
{
	double statistic;
	const Network * pNetwork = this->pNetwork();
	int egoDegree = pNetwork->outDegree(this->ego());
	int alterDegree = pNetwork->inDegree(alter);

	if (this->lroot)
	{
		statistic =
			this->lsqrtTable->sqrt(egoDegree) *
				this->lsqrtTable->sqrt(alterDegree);
	}
	else
	{
		statistic = egoDegree * alterDegree;
	}

	return statistic;
}

}
