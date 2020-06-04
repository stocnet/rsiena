/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveTripletsEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * TransitiveTripletsEffect.
 *****************************************************************************/

#include "TransitiveTripletsEffect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
TransitiveTripletsEffect::TransitiveTripletsEffect(
		const EffectInfo * pEffectInfo, bool twoPath, bool twoInStar) : 
	NetworkEffect(pEffectInfo)
{
	this->ltwoPath = twoPath;
	this->ltwoInStar = twoInStar;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveTripletsEffect::calculateContribution(int alter) const
{
	// If we are introducing a tie from the ego i to the alter j, then each
	// two-path from i to j contributes one transitive triplet, just as each
	// in-star between i and j.
	double contribution = 0;

	if (ltwoPath) contribution += this->pTwoPathTable()->get(alter);
	if (ltwoInStar) contribution += this->pInStarTable()->get(alter);

	return contribution;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double TransitiveTripletsEffect::tieStatistic(int alter)
{
	return this->pTwoPathTable()->get(alter);
}

}
