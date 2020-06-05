/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PrimaryLayer.h
 *
 * Description: This module stores for each actor all its neighbors at
 * distance two with respect to the observed network.
 *****************************************************************************/

#ifndef PrimaryLayer_H_
#define PrimaryLayer_H_

#include "NetworkLayer.h"
#include "../Network.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: PrimaryLayer class
// ----------------------------------------------------------------------------

class PrimaryLayer: public NetworkLayer {
public:
	PrimaryLayer();
	virtual ~PrimaryLayer();
	void onTieIntroductionEvent(const Network& rNetwork, const int ego, const int alter);
	void onTieWithdrawalEvent(const Network& rNetwork, const int ego, const int alter);
	void onNetworkClearEvent(const Network& rNetwork);
	void onNetworkDisposeEvent(const Network& rNetwork);
	int size(int actor); // nl

	const Network * pCounts() const;
	const Network * pLayer() const;

protected:
	void initialize(const Network& rNetwork); // nl

private:
	void initializeOneMode(const Network& rNetwork);
	void initializeTwoMode(const Network& rNetwork);
	void modify2PathCountOneMode(const Network& rNetwork, int ego, int alter, int val);
	void modify2PathCountTwoMode(const Network& rNetwork, int ego, int alter, int val);
	void modifyTieValue(const Network& rNetwork, int ego, int alter, int val);
	void updateSingleTieValue(const Network& rNetwork, int ego, int alter, int val);

	Network* lpCounts;
	Network* lpLayer;

};

} // namespace siena
#endif // PrimaryLayer_H_
