/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DenseTriadsEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * DenseTriadsEffect.
 *****************************************************************************/

#include <stdexcept>
#include "DenseTriadsEffect.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
DenseTriadsEffect::DenseTriadsEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
	this->ldensity = (int) pEffectInfo->internalEffectParameter();

	if (this->ldensity != 5 && this->ldensity != 6)
	{
		throw invalid_argument("Parameter value 5 or 6 expected.");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double DenseTriadsEffect::calculateContribution(int alter) const
{
	int change = 0;

	// It is assumed in the comments below, that there's no tie currently
	// from the ego i to the alter j, i.e. the tie is about to be introduced.

	if (this->ldensity == 6)
	{
		if (this->inTieExists(alter))
		{
			// The dyad (i,j) is becoming complete, so there will be a new
			// full triad for each h with reciprocated ties to both i and j.
			// Each such triad contributes twice to s_i; once for each of the
			// ties (i,j) and (i,h).

			change = 2 * this->pRRTable()->get(alter);
		}
	}
	else
	{
		if (this->inTieExists(alter))
		{
			// Okay, we have a complete dyad (i,j), so each h with at least
			// 3 ties to/from i and j will form a triad with at least 5 ties.
			// Consider a third actor h and all possible combinations of
			// at least 3 ties with i and j:
			//    i->h h->i j->h h->j  # of h   # of contributions
			// A: +    +    +    +     RR       1 (i,h,j) was dense before
			// B: +    +    +    -     RB - RR  2
			// C: +    +    -    +     RF - RR  2
			// D: +    -    +    +     FR - RR  2
			// E: -    +    +    +     BR - RR  1 no tie from i to h
			// The column # of h shows how the number of such actors h may
			// be computed. RB stands for the number of ways we can reach j
			// from i by following one reciprocated tie and one backward tie.
			// RF stands for reciprocated-forward, etc. For example, to
			// compute the number of actors h of type B (second row), we take
			// the number of actors h with a reciprocated tie to i and a (not
			// necessarily reciprocated) tie from j, which is given by RB.
			// Now, since RR is contained in RB, we subtract RR to obtain the
			// number of actors h with a reciprocated tie to i, and a
			// non-reciprocated tie from j. Other rows of the table are treated
			// similarly. The last column shows how a triad of such a type
			// contributes to the statistic s_i. All types contribute at least
			// once because the tie (i,j) is being introduced, but some types
			// of triads contribute twice because the tie (i,h) existed before,
			// but the triad (i,h,j) wasn't dense before introduction of the
			// tie (i,j).

			change =  2 * this->pRFTable()->get(alter) +
				2 * this->pRBTable()->get(alter) +
				2 * this->pFRTable()->get(alter) +
				this->pBRTable()->get(alter) -
				6 * this->pRRTable()->get(alter);
		}
		else
		{
			// The dyad (i,j) is not complete, so we need four more ties for
			// a triad to be dense. The number of actors h with reciprocated
			// ties to both i and j (RR) gives the number of such triads.
			// Note, that each of these triads contributes twice to the
			// statistic s_i, once per each tie (i,j) and (i,h).

			change = 2 * this->pRRTable()->get(alter);
		}
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double DenseTriadsEffect::tieStatistic(int j)
{
	int statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();
	const OneModeNetwork * pOneModeNetwork =
		dynamic_cast<const OneModeNetwork *>(pNetwork);

	// Get the number of dense triads involving the tie (i,j).

	if ((this->ldensity == 6 && this->lmark[j] == this->lbaseMark + 2) ||
		(this->ldensity == 5 && this->lmark[j] == this->lbaseMark + 1))
	{
		// We need complete dyads (j,h) and (i,h) to have
		// a dense triad.

		// Iterate over complete dyads (j,h)

		for (CommonNeighborIterator iterH =
				pOneModeNetwork->reciprocatedTies(j);
			iterH.valid();
			iterH.next())
		{
			int h = iterH.actor();

			// Test if the dyad (i,h) is complete

			if (this->lmark[h] == this->lbaseMark + 2)
			{
				statistic++;
			}
		}
	}
	else if (this->ldensity == 5)
	{
		// The dyad (i,j) is complete, because the previous
		// condition would hold otherwise.

		// Iterate over outgoing and incoming ties of j simultaneously.

		IncidentTieIterator outIter = pNetwork->outTies(j);
		IncidentTieIterator inIter = pNetwork->inTies(j);

		while (outIter.valid() || inIter.valid())
		{
			// Get the current out- or in-neighbor (or n, if we
			// have run out of ties).

			int h1 = n;
			int h2 = n;

			if (outIter.valid())
			{
				h1 = outIter.actor();
			}

			if (inIter.valid())
			{
				h2 = inIter.actor();
			}

			if (h1 == h2)
			{
				// The dyad (j,h1) is complete, so we need just
				// one tie between i and h1.

				if (this->lmark[h1] > this->lbaseMark)
				{
					statistic++;
				}

				outIter.next();
				inIter.next();
			}
			else if (h1 < h2)
			{
				// The dyad (j,h1) has only one tie, so we need
				// a complete dyad (i,h1).

				if (this->lmark[h1] == this->lbaseMark + 2)
				{
					statistic++;
				}

				outIter.next();
			}
			else
			{
				// The dyad (j,h2) has only one tie, so we need
				// a complete dyad (i,h2).

				if (this->lmark[h2] == this->lbaseMark + 2)
				{
					statistic++;
				}

				inIter.next();
			}
		}
	}

	return statistic;
}


/**
 * This method is called at the start of the calculation of the statistic.
 */
void DenseTriadsEffect::initializeStatisticCalculation()
{
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();

	// Allocate the helper array of marks

	this->lmark = new int[n];

	for (int i = 0; i < n; i++)
	{
		this->lmark[i] = 0;
	}
}


/**
 * This method is called at the end of the calculation of the statistic.
 */
void DenseTriadsEffect::cleanupStatisticCalculation()
{
	delete[] this->lmark;
}


/**
 * This method is called right before summing up the contributions of the
 * outgoing ties of the given ego in the calculation of the statistic.
 */
void DenseTriadsEffect::onNextEgo(int i)
{
	const Network * pNetwork = this->pNetwork();
	this->lbaseMark = 2 * i;

	// Count for each actor h the number of ties between i and h
	// (0, 1, or 2). We represent these number as:
	// mark[h] = baseMark + 2 if there are mutual ties between i and h,
	// mark[h] = baseMark + 1 if only one of the mutual ties is present,
	// mark[h] <= baseMark otherwise.

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

}
