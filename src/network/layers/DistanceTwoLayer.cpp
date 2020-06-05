/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DistanceTwoLayer.cpp
 *
 * Description: This module stores for each actor all its neighbors at
 * distance two with respect to the observed network.
 *****************************************************************************/

#include "DistanceTwoLayer.h"

#include <vector>
#include <math.h>

#include "../IncidentTieIterator.h"
#include "../iterators/UnionTieIterator.h"
#include "../iterators/AdvUnionTieIterator.h"

namespace siena {

using namespace std;

typedef map<int, int> TieMap;

/**
 * Constructor.
 */
DistanceTwoLayer::DistanceTwoLayer() :
		NetworkLayer(), //
		lpAdjacencies(0) {
}

/**
 * Destructor.
 */
DistanceTwoLayer::~DistanceTwoLayer() {
	if (lpAdjacencies != 0) {
	delete[] lpAdjacencies;
}
}

/**
 * @copydoc INetworkChangeListener::onTieIntroductionEvent()
 */
void DistanceTwoLayer::onTieIntroductionEvent(const Network& rNetwork,
		const int ego, const int alter) {
	if (rNetwork.isOneMode()) {
		modify2PathCountOneMode(rNetwork, ego, alter, 1);
	} else {
		modify2PathCountTwoMode(rNetwork, ego, alter, -1);
	}
}

/**
 * @copydoc INetworkChangeListener::onTieWithdrawalEvent()
 */
void DistanceTwoLayer::onTieWithdrawalEvent(const Network& rNetwork,
		const int ego, const int alter) {
	if (rNetwork.isOneMode()) {
		modify2PathCountOneMode(rNetwork, ego, alter, -1);
	} else {
		modify2PathCountTwoMode(rNetwork, ego, alter, -1);
	}
}

/**
 * @copydoc INetworkChangeListener::onNetworkClearEvent()
 */
void DistanceTwoLayer::onNetworkClearEvent(const Network& rNetwork) {
	for (int i = 0; i < rNetwork.n(); ++i) {
		lpAdjacencies[i].clear();
	}
}

/**
 * Returns the actor's neighbors at distance two.
 */
IncidentTieIterator DistanceTwoLayer::getDistanceTwoNeighbors(int ego) const {
	return IncidentTieIterator(lpAdjacencies[ego]);
}

/**
 * @copydoc NetworkLayer::initialize()
 */
void DistanceTwoLayer::initialize(const Network& rNetwork) {
	lpAdjacencies = new map<int, int> [rNetwork.n()];
	if (rNetwork.isOneMode()) {
		initializeOneMode(rNetwork);
	} else {
		initializeTwoMode(rNetwork);
	}
}

/**
 * Initializes the layer given the reference network is a one mode
 * network.
 */
void DistanceTwoLayer::initializeOneMode(const Network& rNetwork) {
	for (int i = 0; i < rNetwork.n(); ++i) {
		std::vector<int> neighAtDistOne;
		// avoid the time to copy
		neighAtDistOne.reserve(rNetwork.outDegree(i) + rNetwork.inDegree(i));
		// we could do this all with UnionIterators but it is much slower
		IncidentTieIterator inIter = rNetwork.inTies(i);
		IncidentTieIterator outIter = rNetwork.outTies(i);
		for (UnionTieIterator iter(inIter, outIter); iter.valid();
				iter.next()) {
			// take care of loops
			if (iter.actor() != i) {
				neighAtDistOne.push_back(iter.actor());
			}
		}
		// construct all pairs
		vector<int>::const_iterator iterEnd = neighAtDistOne.end();
		for (vector<int>::const_iterator outerIter = neighAtDistOne.begin();
				outerIter != iterEnd; ++outerIter) {
			int ego = *outerIter;
			for (vector<int>::const_iterator innerIter = outerIter + 1;
					innerIter != iterEnd; ++innerIter) {
				modifyTieValue(ego, *innerIter, 1);
			}
		}
	}
}

/**
 * Initializes the layer given the reference network is a two mode
 * network.
 */
void DistanceTwoLayer::initializeTwoMode(const Network& rNetwork) {
	// this is a two mode network so we do not need to check for loops
	// nor do we have to store the reciever two paths.
	for (int i = 0; i < rNetwork.m(); ++i) {
		// construct all pairs
		for (IncidentTieIterator outerIter = rNetwork.inTies(i);
				outerIter.valid(); outerIter.next()) {
			int outerActor = outerIter.actor();
			// copy the iterator
			IncidentTieIterator innerIter(outerIter);
			// move to the next position
			innerIter.next();
			for (; innerIter.valid(); innerIter.next()) {
				modifyTieValue(outerActor, innerIter.actor(), 1);
			}
		}
	}
}

/**
 * Modifies the two-path count given the case the observed network is a
 * one mode network.
 * @param[in] rNetwork The observed network
 * @param[in] ego The ego of the modified tie
 * @param[in] alter The alter of the modified tie
 * @param[in[ val The magnitude of modification
 */
void DistanceTwoLayer::modify2PathCountOneMode(const Network& rNetwork, int ego,
		int alter, int val) {
	// if it is a loop or the edge (alter,ego) exists we have nothing to do
	if (ego == alter || rNetwork.hasEdge(alter, ego)) {
		return;
	}
	IncidentTieIterator inIter = rNetwork.inTies(ego);
	IncidentTieIterator outIter = rNetwork.outTies(ego);
	UnionTieIterator egoIter(inIter, outIter);
	inIter = rNetwork.inTies(alter);
	outIter = rNetwork.outTies(alter);
	UnionTieIterator alterIter(inIter, outIter);
	AdvUnionTieIterator iter(ego, alter, egoIter, alterIter);
	for (; iter.valid();
			iter.next()) {
		int curNeighbor = iter.actor();
		// check whether the current neighbor is ego or alter itself
		if (curNeighbor != ego && curNeighbor != alter) {
			// if it is a common neighbor it creates a new 2-path with both
			if (iter.isCommon()) {
				modifyTieValue(curNeighbor, ego, val);
				modifyTieValue(curNeighbor, alter, val);
			} else {
				// get the inactive iter id (ego or alter)
				int inactiveIterID = iter.getInactiveIterID();
				modifyTieValue(curNeighbor, inactiveIterID, val);
			}
		}
	}
}

/**
 * Modifies the two-path count given the case the observed network is a
 * two mode network.
 * @param[in] rNetwork The observed network
 * @param[in] ego The ego of the modified tie
 * @param[in] alter The alter of the modified tie
 * @param[in[ val The magnitude of modification
 */
void DistanceTwoLayer::modify2PathCountTwoMode(const Network& rNetwork, int ego,
		int alter, int val) {
	// in a two mode network the exist no triangles, therefore it is
	// sufficient to iterate over all incoming ties of alter
	for (IncidentTieIterator iter = rNetwork.inTies(alter); iter.valid();
			iter.next()) {
		if (iter.actor() != ego) {
			modifyTieValue(ego, iter.actor(), val);
		}
	}
}

/**
 * Updates the tie value of <i>(ego,alter)</i> and <i>(alter,ego)</i>
 * by <i>val</i>
 * @param[in] ego The tie's ego
 * @param[in] alter The tie's alter
 * @param[in] val The magnitude of modification
 */
void DistanceTwoLayer::modifyTieValue(int ego, int alter, int val) {
	updateSingleTieValue(ego, alter, val);
	updateSingleTieValue(alter, ego, val);
}

void DistanceTwoLayer::onNetworkDisposeEvent(const Network& rNetwork) {
	// clear everything
	clear(rNetwork.n());
}

/**
 * Updates the value of tie <i>(ego,alter)</i> by <i>val</i>.
 * @param[in] ego The ego of the modified tie
 * @param[in] alter The alter of the modified tie
 * @param[in[ val The magnitude of modification
 */
void DistanceTwoLayer::updateSingleTieValue(int ego, int alter, int val) {
	TieMap& egoMap = lpAdjacencies[ego];
	TieMap::iterator iter = egoMap.lower_bound(alter);
	// if we found the element update the value and if needed remove
	// the edge from the layer
	if (iter != egoMap.end() && !egoMap.key_comp()(alter, iter->first)) {
		int newVal = iter->second + val;
		if (newVal) {
			iter->second = newVal;
		} else {
			egoMap.erase(iter);
		}
	} else {
		// insert the edge and use iter as hint (saves log)
		egoMap.insert(iter, TieMap::value_type(alter, val));
	}
}

int DistanceTwoLayer::size(int actor) {
	return lpAdjacencies[actor].size();
}

void siena::DistanceTwoLayer::clear(int numOfActors) {
	// TODO: probably not needed
	for (int i = 0; i < numOfActors; i++) {
		lpAdjacencies[i].clear();
	}
	// delete array
	delete[] lpAdjacencies;
	lpAdjacencies = 0;
}

} /* namespace siena */
