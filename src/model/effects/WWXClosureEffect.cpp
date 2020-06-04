/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WWXClosureEffect.cpp
 *
 * Description: This file contains the implementation of the
 * WWXClosureEffect class.
 *****************************************************************************/

#include "WWXClosureEffect.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "network/IncidentTieIterator.h"
#include "data/DyadicCovariateValueIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] out1 and out2 indicate directionality of W ties
 */
WWXClosureEffect::WWXClosureEffect(const EffectInfo * pEffectInfo,
		bool out1, bool out2) :
	DyadicCovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsums = 0;
	this->lout1 = out1;
	this->lout2 = out2;
}


/**
 * Destructor.
 */
WWXClosureEffect::~WWXClosureEffect()
{
	delete[] this->lsums;
	this->lsums = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void WWXClosureEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	DyadicCovariateDependentNetworkEffect::initialize(pData,
		pState,
		period,
		pCache);

	delete[] this->lsums;
	this->lsums = new double[this->pNetwork()->n()];
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void WWXClosureEffect::preprocessEgo(int ego)
{
	DyadicCovariateDependentNetworkEffect::preprocessEgo(ego);
	this->calculateSums(ego, this->pNetwork(), this->lsums);
}


/**
 * For each j and the given i, this method calculates the sum
 * sum_h w_{ih} w_{hj}, or other directionalities depending on out1 and out2.
 */
void WWXClosureEffect::calculateSums(int i,
	const Network * pNetwork,
	double * sums) const
{
	int n = pNetwork->n();

	// Initialize

	for (int j = 0; j < n; j++)
	{
		sums[j] = 0;
	}


	if (this->lout1)
	// Iterate over all h with non-zero non-missing w_{ih}
	{
		for (DyadicCovariateValueIterator iterH = this->rowValues(i);
			iterH.valid();
			iterH.next())
		{
			int h = iterH.actor();
	
			if (this->lout2)
			{

		// Iterate over all j with non-zero non-missing w_{hj}
				for (DyadicCovariateValueIterator iterJ = this->rowValues(h);
					iterJ.valid();
					iterJ.next())
				{
					int j = iterJ.actor();
	
					// Add the term w_{ih} w_{hj}
					sums[j] += iterH.value() * iterJ.value();
				}
			}
			else
			{
		// Iterate over all j with non-zero non-missing w_{jh}
				for (DyadicCovariateValueIterator iterJ = this->columnValues(h);
					iterJ.valid();
					iterJ.next())
				{
					int j = iterJ.actor();
	
					// Add the term w_{ih} w_{jh}
					sums[j] += iterH.value() * iterJ.value();
				}
			}
		}
	}
	else
	// Iterate over all h with non-zero non-missing w_{hi}
	{
		for (DyadicCovariateValueIterator iterH = this->columnValues(i);
			iterH.valid();
			iterH.next())
		{
			int h = iterH.actor();
			if (this->lout2)
			{
		// Iterate over all j with non-zero non-missing w_{hj}
				for (DyadicCovariateValueIterator iterJ = this->rowValues(h);
					iterJ.valid();
					iterJ.next())
				{
					int j = iterJ.actor();
	
					// Add the term w_{hi} w_{hj}
					sums[j] += iterH.value() * iterJ.value();
				}
			}
			else
			{
		// Iterate over all j with non-zero non-missing w_{jh}
				for (DyadicCovariateValueIterator iterJ = this->columnValues(h);
					iterJ.valid();
					iterJ.next())
				{
					int j = iterJ.actor();
	
					// Add the term w_{hi} w_{jh}
					sums[j] += iterH.value() * iterJ.value();
				}
			}
		}
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double WWXClosureEffect::calculateContribution(int alter) const
{
	return this->lsums[alter];
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double WWXClosureEffect::tieStatistic(int alter)
{
	return this->lsums[alter];
}

}
