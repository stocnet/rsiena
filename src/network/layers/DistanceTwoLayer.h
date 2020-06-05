/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DistanceTwoLayer.h
 *
 * Description: This module stores for each actor all its neighbors at
 * distance two with respect to the observed network.
 *****************************************************************************/

#ifndef DISTANCETWOLAYER_H_
#define DISTANCETWOLAYER_H_

#include <map>

#include "NetworkLayer.h"
#include "../Network.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: DistanceTwoLayer class
// ----------------------------------------------------------------------------

class DistanceTwoLayer: public NetworkLayer {
public:

	DistanceTwoLayer();

	virtual ~DistanceTwoLayer();

	void onTieIntroductionEvent(const Network& rNetwork, const int ego,
			const int alter);

	void onTieWithdrawalEvent(const Network& rNetwork, const int ego,
			const int alter);

	void onNetworkClearEvent(const Network& rNetwork);

	void onNetworkDisposeEvent(const Network& rNetwork);

	IncidentTieIterator getDistanceTwoNeighbors(int ego) const;

	virtual int size(int actor);

	void clear(int numOfActors);

protected:

	void initialize(const Network& rNetwork);

private:

	void initializeOneMode(const Network& rNetwork);

	void initializeTwoMode(const Network& rNetwork);

	void modify2PathCountOneMode(const Network& rNetwork, int ego, int alter,
			int val);

	void modify2PathCountTwoMode(const Network& rNetwork, int ego, int alter,
			int val);

	void modifyTieValue(int ego, int alter, int val);

	void updateSingleTieValue(int ego, int alter, int val);

	std::map<int, int>* lpAdjacencies;
};

} /* namespace siena */
#endif /* DISTANCETWOLAYER_H_ */
