/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PrimaryLayer.cpp
 *
 * Description: This module stores for each actor all its neighbors at
 * distance two with respect to the observed network.
 *****************************************************************************/

 
#include <R_ext/Error.h>
#include "PrimaryLayer.h"

#include <vector>
#include <math.h>

#include "network/OneModeNetwork.h"

#include "../IncidentTieIterator.h"
#include "../iterators/UnionTieIterator.h"
#include "../iterators/AdvUnionTieIterator.h"
#include <Rinternals.h>


namespace siena {

using namespace std;

// typedef map<int, int> TieMap;

PrimaryLayer::PrimaryLayer() :
		NetworkLayer(), //
		lpCounts(0),
		lpLayer(0) {
}

PrimaryLayer::~PrimaryLayer() {
	delete lpCounts;
	delete lpLayer;
	lpCounts = 0;
	lpLayer = 0;
}

void PrimaryLayer::onTieIntroductionEvent(const Network& rNetwork, const int ego, const int alter) {
	// LOGS(Priority::ERROR) << "add: " << ego << "," << alter;
	if (rNetwork.isOneMode()) {
		modify2PathCountOneMode(rNetwork, ego, alter, 1);
	} else {
		modify2PathCountTwoMode(rNetwork, ego, alter, -1);
	}
}

/**
 * @copydoc INetworkChangeListener::onTieWithdrawalEvent()
 */
void PrimaryLayer::onTieWithdrawalEvent(const Network& rNetwork, const int ego, const int alter) {
	// LOGS(Priority::ERROR) << "del: " << ego << "," << alter;
	if (rNetwork.isOneMode()) {
		modify2PathCountOneMode(rNetwork, ego, alter, -1);
	} else {
		modify2PathCountTwoMode(rNetwork, ego, alter, -1);
	}
}

/**
 * @copydoc INetworkChangeListener::onNetworkClearEvent()
 */
void PrimaryLayer::onNetworkClearEvent(const Network& /*rNetwork*/) {
	lpCounts->clear();
	lpLayer->clear();
}

/**
 * @copydoc NetworkLayer::initialize()
 */
void PrimaryLayer::initialize(const Network& rNetwork) {
	// LOGS(Priority::ERROR) << "init";
	// We should not dispose here!  The network is hold by state; it is set only
	// in its ctor but initialized/disposed multiple times.
	if (lpLayer == 0) {
		if (rNetwork.isOneMode()) {
			lpCounts = new OneModeNetwork(rNetwork.n());
			lpLayer = new OneModeNetwork(rNetwork.n());
		} else {
			Rf_error("not implemented");
		}
	}

	lpCounts->clear();
	lpLayer->clear();

	if (rNetwork.isOneMode()) {
		initializeOneMode(rNetwork);
	} else {
		initializeTwoMode(rNetwork);
	}
	// LOGS(Priority::ERROR) << "init done";
}

/**
 * Initializes the layer given the reference network is a one mode
 * network.
 */
void PrimaryLayer::initializeOneMode(const Network& rNetwork) {
	typedef vector<int>::const_iterator viItr;

	for (int i = 0; i < rNetwork.n(); ++i) {
		std::vector<int> neighAtDistOne;

		// avoid the time to copy
		neighAtDistOne.reserve(rNetwork.outDegree(i) + rNetwork.inDegree(i));

		// collect the undirected neighborhood
		IncidentTieIterator inIter = rNetwork.inTies(i);
		IncidentTieIterator outIter = rNetwork.outTies(i);
		for (UnionTieIterator iter(inIter, outIter); iter.valid(); iter.next()) {
			if (iter.actor() == i) continue; // skip loops

			lpLayer->setTieValue(i, iter.actor(), 1); // copy original network to layer

			// LOGS(Priority::ERROR) << "i " << i << " " << iter.actor();
			// modifyTieValue(rNetwork, i, iter.actor(), 1);
			neighAtDistOne.push_back(iter.actor());
		}

		// construct all pairs
		viItr iterEnd = neighAtDistOne.end();
		for (viItr outerIter = neighAtDistOne.begin(); outerIter != iterEnd; ++outerIter) {
			int ego = *outerIter;
			for (viItr innerIter = outerIter + 1; innerIter != iterEnd; ++innerIter) {
				// LOGS(Priority::ERROR) << "ii " << ego << " " << *innerIter;
				modifyTieValue(rNetwork, ego, *innerIter, 1);
			}
		}
	}
}

/**
 * Initializes the layer given the reference network is a two mode
 * network.
 */
void PrimaryLayer::initializeTwoMode(const Network& rNetwork) {
	Rf_error("primary layer not implemented for two-mode");
}

void PrimaryLayer::modify2PathCountOneMode(const Network& rNetwork, int ego, int alter, int val) {
	if (ego == alter || rNetwork.hasEdge(alter, ego)) return;
	// LOGS(Priority::ERROR) << "init mod_" << ego << "," << alter << " + " << val;

	// undirected neighborhood of ego
	IncidentTieIterator inIter = rNetwork.inTies(ego);
	IncidentTieIterator outIter = rNetwork.outTies(ego);
	UnionTieIterator egoIter(inIter, outIter);
	// undirected neighborhood of alter
	inIter = rNetwork.inTies(alter);
	outIter = rNetwork.outTies(alter);
	UnionTieIterator alterIter(inIter, outIter);
	// union tracking commonness with regard to ego/alter
	AdvUnionTieIterator iter(ego, alter, egoIter, alterIter);

	// for all undirected neightbors of either ends
	for (; iter.valid(); iter.next()) {
		int curNeighbor = iter.actor();
		// check whether the current neighbor is ego or alter itself
		if (curNeighbor == ego || curNeighbor == alter) continue;

		// if it is a common neighbor it creates a new 2-path with both
		if (iter.isCommon()) {
			modifyTieValue(rNetwork, curNeighbor, ego, val);
			modifyTieValue(rNetwork, curNeighbor, alter, val);
		} else {
			int inactiveIterId = iter.getInactiveIterID();
			// get the inactive iter id (ego or alter)
			modifyTieValue(rNetwork, curNeighbor, inactiveIterId, val);
		}
	}

	// might have removed tie ego-alter without changeing 2path counts, we need
	// to check wether the tie is still in the network (had any 2paths)
	modifyTieValue(rNetwork, ego, alter, 0);
}

/**
 * Modifies the two-path count given the case the observed network is a
 * two mode network.
 *
 * @param[in] rNetwork The observed network
 * @param[in] ego The ego of the modified tie
 * @param[in] alter The alter of the modified tie
 * @param[in[ val The magnitude of modification
 */
void PrimaryLayer::modify2PathCountTwoMode(const Network& rNetwork, int ego, int alter, int val) {
	Rf_error("not implemented");
}

void PrimaryLayer::onNetworkDisposeEvent(const Network& /*rNetwork*/) {
	// We should not dispose here!  The network is hold by state; it is set only
	// in its ctor but initialized/disposed multiple times.
	//delete lpLayer;
	//lpLayer = 0;
}

/**
 * Updates the tie value of <i>(ego,alter)</i> and <i>(alter,ego)</i>
 * by <i>val</i>
 *
 * @param[in] ego The tie's ego
 * @param[in] alter The tie's alter
 * @param[in] val The magnitude of modification
 */
void PrimaryLayer::modifyTieValue(const Network& rNetwork, int ego, int alter, int val) {
	updateSingleTieValue(rNetwork, ego, alter, val);
	updateSingleTieValue(rNetwork, alter, ego, val);
}

/**
 * Updates the value of tie <i>(ego,alter)</i> by <i>val</i>.
 * Only updates the store, no logic involved.
 *
 * @param[in] ego The ego of the modified tie
 * @param[in] alter The alter of the modified tie
 * @param[in[ val The magnitude of modification
 */
void PrimaryLayer::updateSingleTieValue(const Network& rNetwork, int ego, int alter, int val) {
	const int cnt = lpCounts->tieValue(ego, alter);
	const int newCnt = cnt + val;
	lpCounts->setTieValue(ego, alter, newCnt);
	if (newCnt > 0 ||
			rNetwork.tieValue(ego, alter) + rNetwork.tieValue(alter, ego) > 0) {
		// LOGS(Priority::ERROR) << " ++ " << ego << "," << alter << " (" << cnt << " -> " << newCnt << ")";
		lpLayer->setTieValue(ego, alter, 1);
	} else {
		// LOGS(Priority::ERROR) << " -- " << ego << "," << alter;
		lpLayer->setTieValue(ego, alter, 0);
	}
}

const Network * PrimaryLayer::pCounts() const {
	return this->lpCounts;
}

const Network * PrimaryLayer::pLayer() const {
	return this->lpLayer;
}

int PrimaryLayer::size(int actor) {
	return lpLayer->outDegree(actor);
}

} // namespace siena
