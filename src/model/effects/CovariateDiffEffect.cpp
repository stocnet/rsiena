/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDiffEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDiffEffect class.
 *****************************************************************************/

#include <cmath>
#include "CovariateDiffEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
CovariateDiffEffect::CovariateDiffEffect(const EffectInfo * pEffectInfo,
		bool diff, int trafo) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
	this->ldiff = diff;
	this->lsquared = (trafo == 2);
	this->labs = (trafo == 1);
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateDiffEffect::calculateContribution(int alter) const
{
	double change = 0;

	if (this->ldiff)
	{
		change = this->value(alter) - this->value(this->ego());
		if (this->lsquared)
		{
			change *= change;
		}
		if (this->labs)
		{
			change = fabs(change);
		}
	}
	else
	{
		if (this->lsquared)
		{
			change = (this->value(alter) * this->value(alter)) + 
						(this->value(this->ego()) * this->value(this->ego()));
		}
		else
		{
			change = this->value(alter) + this->value(this->ego());
		}
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateDiffEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!(this->missing(alter) || this->missing(this->ego())))
	{
		if (this->ldiff)
		{
			statistic = this->value(alter) - this->value(this->ego());
			if (this->labs)
			{
				statistic = fabs(statistic);
			}

			if (this->lsquared)
			{
				statistic *= statistic;
			}
		}
		else
		{
			if (this->lsquared)
			{
				statistic = (this->value(alter) * this->value(alter)) + 
								(this->value(this->ego()) * this->value(this->ego()));
			}
			else
			{
				statistic = this->value(alter) + this->value(this->ego());
			}
		}
	}

	return statistic;
}

}
