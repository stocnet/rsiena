/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Network.h
 *
 * Description: This module defines the class Network for storing directed
 * valued networks.
 *****************************************************************************/

#ifndef NETWORK_H_
#define NETWORK_H_

#include <map>
#include <list>
#include <string>

namespace siena {

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class TieIterator;
class IncidentTieIterator;
class INetworkChangeListener;

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

enum ChangeType {
	REPLACE, INCREASE
};

// ----------------------------------------------------------------------------
// Section: Network class
// ----------------------------------------------------------------------------

/**
 * This class defines a directed network on two sets of actors, namely, those
 * acting as senders and receivers of ties, respectively. A single set of
 * actors may act as both senders and receivers, in which case the use of
 * OneModeNetwork is recommended. The ties are valued and multiple ties between
 * the same pair of actors are forbiden.
 */
class Network {
public:
	Network(int n, int m);
	Network(const Network & rNetwork);
	Network & operator=(const Network & rNetwork);
	virtual Network * clone() const;
	virtual ~Network();

	int n() const;
	int m() const;
	int tieCount() const;

	void setTieValue(int i, int j, int v);
	int tieValue(int i, int j) const;
	int increaseTieValue(int i, int j, int v);
	virtual void clear();
	void clearInTies(int actor);
	void clearOutTies(int actor);

	TieIterator ties() const;
	IncidentTieIterator inTies(int i) const;
	IncidentTieIterator outTies(int i) const;
	IncidentTieIterator inTies(int i, std::string mess) const;
	IncidentTieIterator inTies(int i, int lowerBound) const;
	IncidentTieIterator outTies(int i, int lowerBound) const;

	int inDegree(int i) const;
	int outDegree(int i) const;

	int minTieValue() const;
	int maxTieValue() const;

	bool complete() const;
	bool hasEdge(int ego, int alter) const;
	virtual bool isOneMode() const;

	int outTwoStarCount(int i, int j) const;
	int inTwoStarCount(int i, int j) const;

	void addNetworkChangeListener(INetworkChangeListener* const listener);

	void removeNetworkChangeListener(INetworkChangeListener* const listener);

	inline int modificationCount() const;

protected:
	virtual int changeTieValue(int i, int j, int v, ChangeType type);
	virtual void onTieWithdrawal(int i, int j);
	virtual void onTieIntroduction(int i, int j);
	void checkSenderRange(int i) const;
	void checkReceiverRange(int i) const;
	void checkReceiverRange(int i, std::string message) const;
	virtual int maxTieCount() const;

	// set of network change listener
	std::list<INetworkChangeListener*> lNetworkChangeListener;
private:
	void allocateArrays();
	void deleteArrays();
	void fireNetworkDisposeEvent();
	void fireNetworkClearEvent() const;
	void fireIntroductionEvent(int ego, int alter) const;
	void fireWithdrawalEvent(int ego, int alter) const;

	// The number of senders
	int ln {};

	// The number of receivers
	int lm {};

	// An array of maps storing outgoing ties of each sender. A tie (i,j)
	// with a non-zero value v is stored as a pair (j,v) in lpOutTies[i].

	std::map<int, int> * lpOutTies;

	// An array of maps storing incoming ties of each receiver. A tie (i,j)
	// with a non-zero value v is stored as a pair (i,v) in lpInTies[j].

	std::map<int, int> * lpInTies;

	// The number of ties of this network
	int ltieCount {};

	// This variable is initially 0 and incremented each time the network
	// is changed.

	int lmodificationCount {};

};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the number of times this network has been changed.
 */
int Network::modificationCount() const {
	return this->lmodificationCount;
}

}

#endif /*NETWORK_H_*/
