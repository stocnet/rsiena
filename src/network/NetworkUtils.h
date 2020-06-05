/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkUtils.h
 *
 * Description: This module contains some utilities specific to the
 * 'network' library.
 *****************************************************************************/

#include <vector>

#ifndef NETWORKUTILS_H_
#define NETWORKUTILS_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

/**
 * This enumeration defines two possible directions of traversing a tie.
 */
enum Direction {FORWARD, BACKWARD, RECIPROCAL};


// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class IncidentTieIterator;
class Network;


// ----------------------------------------------------------------------------
// Section: Utility functions
// ----------------------------------------------------------------------------

int commonActorCount(IncidentTieIterator iter1, IncidentTieIterator iter2);
Network * symmetricDifference(const Network * pNetwork1,
	const Network * pNetwork2);

void subtractNetwork(Network * pNetwork,
	 const Network * pMissingTieNetwork);

void replaceNetwork(Network * pNetwork,
	const Network * pValueNetwork, const Network * pDecisionNetwork);

std::vector<int> * primarySetting(const Network * pNetwork, int ego);
}

#endif /*NETWORKUTILS_H_*/
