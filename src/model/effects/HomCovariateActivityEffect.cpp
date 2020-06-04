/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HomCovariateActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * HomCovariateActivityEffect class.
 *****************************************************************************/

#include <cmath>
#include "HomCovariateActivityEffect.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
HomCovariateActivityEffect::HomCovariateActivityEffect(const EffectInfo * pEffectInfo,
	bool same) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsame = same;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double HomCovariateActivityEffect::calculateContribution(int alter) const
{
	double altervalue = this->value(alter);
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();
	if (lsame)
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
		{
			// Get the receiver of the outgoing tie.
			int h = iter.actor();
			if ((h != alter) && (fabs(this->value(h) - altervalue) < EPSILON))
			{
				contribution++;
			}
		}
	}
	else
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
		{
			// Get the receiver of the outgoing tie.
			int h = iter.actor();
			if ((h != alter) && (fabs(this->value(h) - altervalue) >= EPSILON))
			{
				contribution++;
			}
		}
	}
	contribution *= 2;
	contribution++;
	return contribution;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double HomCovariateActivityEffect::tieStatistic(int alter)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (!(this->missing(alter)))
	{
		double altervalue = this->value(alter);
		if (lsame)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
				iter.valid();
				iter.next())
			{
				// Get the receiver of the outgoing tie.
				int h = iter.actor();
				if ((!this->missing(h)) &&
						(fabs(this->value(h) - altervalue) < EPSILON))
				{
					contribution++;
				}
			}
		}
		else
		{
			for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
				iter.valid();
				iter.next())
			{
				// Get the receiver of the outgoing tie.
				int h = iter.actor();
				if ((!this->missing(h)) &&
						(fabs(this->value(h) - altervalue) >= EPSILON))
				{
					contribution++;
				}
			}
		}
	}
	return contribution;
}

}
