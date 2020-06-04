/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ThreeCyclesEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * ThreeCyclesEffect.
 *****************************************************************************/

#include "ThreeCyclesEffect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
ThreeCyclesEffect::ThreeCyclesEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double ThreeCyclesEffect::calculateContribution(int alter) const
{
	// The absolute value of the contribution is going to be the number
	// of two-paths from the alter to the ego -- a number, which is stored
	// in the table of reverse two paths.
	return this->pReverseTwoPathTable()->get(alter);
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double ThreeCyclesEffect::tieStatistic(int alter)
{
	return this->pReverseTwoPathTable()->get(alter) / 3.0;
}

}
