/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkLayer.h
 *
 * Description: This module defines the abstract NetworkLayer implementing
 * the INetworkChangeListener class. It serves as the base class for any
 * network layer.
 *****************************************************************************/

#ifndef NETWORKLAYER_H_
#define NETWORKLAYER_H_

#include "../INetworkChangeListener.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: NetworkLayer abstract class
// ----------------------------------------------------------------------------

class NetworkLayer: public INetworkChangeListener {
public:

	/**
	 * Destructor.
	 */
	virtual ~NetworkLayer() {
	}

	void onInitializationEvent(const Network& rNetwork) {
		initialize(rNetwork);
	}

	virtual int size(int actor) = 0;

protected:

	/**
	 * Constructor.
	 */
	NetworkLayer() :
			INetworkChangeListener() {
	}

	/**
	 * Initializes the network layer.
	 */
	virtual void initialize(const Network& rNetwork) = 0;

private:
	// disable copy constructor and copy assignment
	NetworkLayer& operator=(const NetworkLayer& rhs);
	NetworkLayer(const NetworkLayer& rhs);
};

} /* namespace siena */

#endif /* NETWORKLAYER_H_ */
