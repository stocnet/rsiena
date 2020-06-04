/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DenseTriadsSimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DenseTriadsSimilarityEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <cstdlib>
#include "DenseTriadsSimilarityEffect.h"
#include "model/EffectInfo.h"
#include "network/OneModeNetwork.h"
#include "network/CommonNeighborIterator.h"
#include "network/IncidentTieIterator.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
DenseTriadsSimilarityEffect::DenseTriadsSimilarityEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->ldensity = (int) pEffectInfo->internalEffectParameter();
	this->lmark = 0;
	this->lbaseMark = 0;

	if (this->ldensity != 5 && this->ldensity != 6)
	{
		throw invalid_argument("Parameter value 5 or 6 expected.");
	}
}


/**
 * Destructor.
 */
DenseTriadsSimilarityEffect::~DenseTriadsSimilarityEffect()
{
	delete[] this->lmark;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DenseTriadsSimilarityEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkDependentBehaviorEffect::initialize(pData, pState, period, pCache);
	delete[] this->lmark;
	int n = this->pNetwork()->n();
	this->lmark = new int[n];
	this->lbaseMark = 0;

	for (int i = 0; i < n; i++)
	{
		this->lmark[i] = 0;
	}
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double DenseTriadsSimilarityEffect::calculateChangeContribution(int i,
	int difference)
{
	this->updateMarks(i);

	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(
			"One-mode network expected in DenseTriadsBehaviorEffect");
	}

	double totalSimilarityChange = 0;
	int oldValue = this->value(i);
	int newValue = oldValue + difference;

	if (this->ldensity == 6)
	{
		for (CommonNeighborIterator iterI = pNetwork->reciprocatedTies(i);
			iterI.valid();
			iterI.next())
		{
			int j = iterI.actor();

			for (CommonNeighborIterator iterJ = pNetwork->reciprocatedTies(j);
				iterJ.valid();
				iterJ.next())
			{
				int h = iterJ.actor();

				if (this->lmark[h] == this->lbaseMark + 2)
				{
					int alterValue = this->value(j);
					totalSimilarityChange +=
						abs(oldValue - alterValue) -
							abs(newValue - alterValue);

					alterValue = this->value(h);
					totalSimilarityChange +=
						abs(oldValue - alterValue) -
							abs(newValue - alterValue);
				}
			}
		}

		// Each triad was visited twice
		totalSimilarityChange /= 2;
	}
	else
	{
		for (CommonNeighborIterator iterI = pNetwork->reciprocatedTies(i);
			iterI.valid();
			iterI.next())
		{
			int j = iterI.actor();
			IncidentTieIterator outIter = pNetwork->outTies(j);
			IncidentTieIterator inIter = pNetwork->inTies(j);

			while (outIter.valid())
			{
				while (inIter.valid() && inIter.actor() < outIter.actor())
				{
					inIter.next();
				}

				int h = outIter.actor();
				bool reciprocatedJH = inIter.valid() && inIter.actor() == h;
				int tieCountBetweenIH = this->lmark[h] - this->lbaseMark;

				if (reciprocatedJH)
				{
					// One more tie needed.

					if (tieCountBetweenIH > 0)
					{
						// Got the required tie. Make sure that the triad
						// is included only once. If there's a reciprocated tie
						// between i and h, then we'll hit this triad twice,
						// but we count it only when j is the smaller of the
						// two alters.

						if (tieCountBetweenIH != 2 || j < h)
						{
							int alterValue = this->value(j);
							totalSimilarityChange +=
								abs(oldValue - alterValue) -
									abs(newValue - alterValue);

							alterValue = this->value(h);
							totalSimilarityChange +=
								abs(oldValue - alterValue) -
									abs(newValue - alterValue);
						}
					}
				}
				else
				{
					// Two more ties needed.

					if (tieCountBetweenIH == 2)
					{
						int alterValue = this->value(j);
						totalSimilarityChange +=
							abs(oldValue - alterValue) -
								abs(newValue - alterValue);

						alterValue = this->value(h);
						totalSimilarityChange +=
							abs(oldValue - alterValue) -
								abs(newValue - alterValue);
					}
				}

				outIter.next();
			}
		}
	}

	// Remember that the similarities are normalized by dividing with the
	// total range of the variable.

	return totalSimilarityChange / this->range();
}


/**
 * Counts for each actor h the number of ties between i and h
 * (0, 1, or 2). We represent these number as:
 * mark[h] = baseMark + 2 if there are mutual ties between i and h,
 * mark[h] = baseMark + 1 if only one of the mutual ties is present,
 * mark[h] <= baseMark otherwise.
 */
void DenseTriadsSimilarityEffect::updateMarks(int i)
{
	const Network * pNetwork = this->pNetwork();
	this->lbaseMark += 2;

	for (IncidentTieIterator iter = pNetwork->inTies(i);
		iter.valid();
		iter.next())
	{
		this->lmark[iter.actor()] = this->lbaseMark + 1;
	}

	for (IncidentTieIterator iter = pNetwork->outTies(i);
		iter.valid();
		iter.next())
	{
		if (this->lmark[iter.actor()] <= this->lbaseMark)
		{
			this->lmark[iter.actor()] = this->lbaseMark + 1;
		}
		else
		{
			this->lmark[iter.actor()]++;
		}
	}
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double DenseTriadsSimilarityEffect::egoStatistic(int i, double * currentValues)
{
	this->updateMarks(i);

	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pNetwork)
	{
		throw runtime_error(
			"One-mode network expected in DenseTriadsBehaviorEffect");
	}

	double totalSimilarity = 0;

	if (this->ldensity == 6)
	{
		for (CommonNeighborIterator iterI = pNetwork->reciprocatedTies(i);
			iterI.valid();
			iterI.next())
		{
			int j = iterI.actor();

			if (!this->missing(this->period(), j) &&
				!this->missing(this->period() + 1, j))
			{
				for (CommonNeighborIterator iterJ =
						pNetwork->reciprocatedTies(j);
					iterJ.valid();
					iterJ.next())
				{
					int h = iterJ.actor();

					if (this->lmark[h] == this->lbaseMark + 2 &&
						!this->missing(this->period(), h) &&
						!this->missing(this->period() + 1, h))
					{
						totalSimilarity +=
							this->similarity(currentValues[i],
								currentValues[j]) +
							this->similarity(currentValues[i],
								currentValues[h]);
					}
				}
			}
		}

		// Each triad was visited twice
		totalSimilarity /= 2;
	}
	else
	{
		for (CommonNeighborIterator iterI = pNetwork->reciprocatedTies(i);
			iterI.valid();
			iterI.next())
		{
			int j = iterI.actor();

			if (!this->missing(this->period(), j) &&
				!this->missing(this->period() + 1, j))
			{
				IncidentTieIterator outIter = pNetwork->outTies(j);
				IncidentTieIterator inIter = pNetwork->inTies(j);

				while (outIter.valid())
				{
					while (inIter.valid() && inIter.actor() < outIter.actor())
					{
						inIter.next();
					}

					int h = outIter.actor();

					if (!this->missing(this->period(), h) &&
						!this->missing(this->period() + 1, h))
					{
						bool reciprocatedJH =
							inIter.valid() && inIter.actor() == h;
						int tieCountBetweenIH =
							this->lmark[h] - this->lbaseMark;

						if (reciprocatedJH)
						{
							// One more tie needed.

							if (tieCountBetweenIH > 0)
							{
								// Got the required tie. Make sure that the
								// triad is included only once. If there's a
								// reciprocated tie between i and h, then we'll
								// hit this triad twice, but we include it only
								// when j is the smaller of the two alters.

								if (tieCountBetweenIH != 2 || j < h)
								{
									totalSimilarity +=
										this->similarity(currentValues[i],
											currentValues[j]) +
										this->similarity(currentValues[i],
											currentValues[h]);
								}
							}
						}
						else
						{
							// Two more ties needed.

							if (tieCountBetweenIH == 2)
							{
								totalSimilarity +=
									this->similarity(currentValues[i],
										currentValues[j]) +
									this->similarity(currentValues[i],
										currentValues[h]);
							}
						}
					}

					outIter.next();
				}
			}
		}
	}

	return totalSimilarity;
}

}
