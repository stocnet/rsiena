/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: RecipdegreeActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * RecipdegreeActivityEffect class.
 *****************************************************************************/

#include <cmath>
#include <string>
#include <stdexcept>
#include "RecipdegreeActivityEffect.h"
#include "utils/SqrtTable.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 */
RecipdegreeActivityEffect::RecipdegreeActivityEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lroot = (std::abs(pEffectInfo->internalEffectParameter()-2) < 0.001);
	this->lsqrtTable = SqrtTable::instance();
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double RecipdegreeActivityEffect::calculateContribution(int alter) const
{
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pONetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreeBehaviorEffect");
	}


	double change = 0;
	double rdegree = pONetwork->reciprocalDegree(this->ego());

	if (this->lroot)
	{
		if (this->outTieExists(alter))
		{
			rdegree --;
		}
		double difrd = 0;
		double rrdegree = this->lsqrtTable->sqrt(rdegree);
		if (this->inTieExists(alter))
		{
			int outd = this->pNetwork()->outDegree(this->ego());
			if (!(this->outTieExists(alter)))
			{
				outd ++;
			}
			difrd = outd * (this->lsqrtTable->sqrt(rdegree+1) - rrdegree);
		}
		change = rrdegree + difrd;
	}
	else
	{
		if (this->inTieExists(alter))
		{
			rdegree += this->pNetwork()->outDegree(this->ego());
		if (this->outTieExists(alter))
			{
				rdegree --;
			}
			else
			{
				rdegree ++;
			}
	}
		change = rdegree;
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double RecipdegreeActivityEffect::tieStatistic(int alter)
{
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pONetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreeBehaviorEffect");
	}
	if (this->lroot)
	{
		return this->lsqrtTable->sqrt(pONetwork->reciprocalDegree(this->ego()));
	}
	else
	{
		return pONetwork->reciprocalDegree(this->ego());
	}
}

}
