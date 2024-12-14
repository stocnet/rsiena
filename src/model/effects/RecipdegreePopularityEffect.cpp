/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: RecipdegreePopularityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * RecipdegreePopularityEffect class.
 *****************************************************************************/

#include <string>
#include <stdexcept>
#include "RecipdegreePopularityEffect.h"
#include "utils/SqrtTable.h"
//#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
RecipdegreePopularityEffect::RecipdegreePopularityEffect(
	const EffectInfo * pEffectInfo, bool root) : NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double RecipdegreePopularityEffect::calculateContribution(int alter) const
{
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pONetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreePopularityEffect");
	}
	
	int degree = pONetwork->reciprocalDegree(alter);
//	int degree = this->pNetwork()->reciprocalDegree(alter);

	if (this->inTieExists(alter))
	{
		// The reciprocal degree will increase after introducing the tie
		degree++;
	}

	double change = degree;

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(degree);
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double RecipdegreePopularityEffect::tieStatistic(int alter)
{
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pONetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreePopularityEffect");
	}
	
	int degree = pONetwork->reciprocalDegree(alter);
	double statistic;

	if (this->lroot)
	{
		statistic = this->lsqrtTable->sqrt(degree);
	}
	else
	{
		statistic = degree;
	}
	return statistic;
}

}
