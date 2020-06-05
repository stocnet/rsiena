/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: INetworkChangeListener.h
 *
 * Description: This module defines the interface INetworkChangeListener.
 * Any class implementing this interface can be added to a network and gets
 * informed about network changes (tie introduction/withdrawal or clearing).
 *****************************************************************************/

#ifndef INETWORKCHANGELISTENER_H_
#define INETWORKCHANGELISTENER_H_

#include "Network.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: INetworkChangeListener interface
// ----------------------------------------------------------------------------

class INetworkChangeListener {
public:

	/**
	 * Destructor.
	 */
	virtual ~INetworkChangeListener() {
	}

	virtual void onInitializationEvent(const Network& rNetwork) = 0;

	/**
	 * Invoked when an tie is introduced to the network.
	 */
	virtual void onTieIntroductionEvent(const Network& rNetwork, const int ego,
			const int alter) = 0;

	/**
	 * Invoked when an tie is withdrawn from the network.
	 */
	virtual void onTieWithdrawalEvent(const Network& rNetwork, const int ego,
			const int alter) = 0;

	/**
	 * Invoked when the network is cleared.
	 */
	virtual void onNetworkClearEvent(const Network& rNetwork) = 0;

	/**
	 * Invoked when the network is disposed.
	 */
	virtual void onNetworkDisposeEvent(const Network& rNetwork) = 0;

protected:

	/**
	 * Constructor.
	 */
	INetworkChangeListener() {
	}
private:

	// disable copy constructor and copy assignment
	INetworkChangeListener(const INetworkChangeListener& rhs);
	INetworkChangeListener& operator=(const INetworkChangeListener& rhs);
};

}
#endif /* INETWORKCHANGELISTENER_H_ */
