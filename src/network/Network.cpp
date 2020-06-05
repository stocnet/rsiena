/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Network.cpp
 *
 * Description: This module implements the class Network for storing directed
 * valued networks.
 *****************************************************************************/

#include "Network.h"

#include <map>
#include <stdexcept>
#include <limits>
#include <algorithm>

#include "TieIterator.h"
#include "IncidentTieIterator.h"
#include "NetworkUtils.h"
#include "INetworkChangeListener.h"
#include "../utils/Utils.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: Construction and destruction
// ----------------------------------------------------------------------------

/**
 * Creates an empty network with <i>n</i> senders and <i>m</i> receivers.
 */
Network::Network(int n, int m) {
	if (n < 0) {
		throw std::invalid_argument("Negative number of senders specified");
	}

	if (m < 0) {
		throw std::invalid_argument("Negative number of receivers specified");
	}

	// Store the parameters

	this->ln = n;
	this->lm = m;

	// Allocate data structures
	this->allocateArrays();

	// No ties for now
	this->ltieCount = 0;
	this->lmodificationCount = 0;
}

/**
 * Creates a copy of the given network.
 */
Network::Network(const Network & rNetwork) {
	this->ln = rNetwork.ln;
	this->lm = rNetwork.lm;

	// Allocate data structures
	this->allocateArrays();

	// Copy everything from rNetwork

	for (int i = 0; i < this->ln; i++) {
		this->lpOutTies[i].insert(rNetwork.lpOutTies[i].begin(),
				rNetwork.lpOutTies[i].end());
	}

	for (int i = 0; i < this->lm; i++) {
		this->lpInTies[i].insert(rNetwork.lpInTies[i].begin(),
				rNetwork.lpInTies[i].end());
	}

	this->ltieCount = rNetwork.ltieCount;
	this->lmodificationCount = 0;
}

/**
 * Assigns the contents of the given network to this network.
 */
Network & Network::operator=(const Network & rNetwork) {
	if (this != &rNetwork) {
		fireNetworkDisposeEvent();
		// Empty the current network structure.

		for (int i = 0; i < this->ln; i++) {
			this->lpOutTies[i].clear();
		}

		for (int i = 0; i < this->lm; i++) {
			this->lpInTies[i].clear();
		}

		// Store the new size of actor sets

		this->ln = rNetwork.ln;
		this->lm = rNetwork.lm;

		// Reallocate the data structures

		this->deleteArrays();
		this->allocateArrays();

		// Copy everything from rNetwork

		for (int i = 0; i < this->ln; i++) {
			this->lpOutTies[i].insert(rNetwork.lpOutTies[i].begin(),
					rNetwork.lpOutTies[i].end());
		}

		for (int i = 0; i < this->lm; i++) {
			this->lpInTies[i].insert(rNetwork.lpInTies[i].begin(),
					rNetwork.lpInTies[i].end());
		}

		this->ltieCount = rNetwork.ltieCount;
		this->lmodificationCount++;
	}

	if (!isOneMode()) {
		for (std::list<INetworkChangeListener*>::const_iterator iter =
				lNetworkChangeListener.begin();
				iter != lNetworkChangeListener.end(); ++iter) {
			(*iter)->onInitializationEvent(*this);
		}
	}
	return *this;
}

/**
 * This method creates a copy of this network.
 */
Network * Network::clone() const {
	return new Network(*this);
}

/**
 * This method allocates the memory for maps of incident ties and
 * various arrays of counters.
 */
void Network::allocateArrays() {
	// Allocate data structures

	this->lpOutTies = new std::map<int, int>[this->ln];
	this->lpInTies = new std::map<int, int>[this->lm];
}

/**
 * Deallocates various arrays used by this network.
 */
void Network::deleteArrays() {
	delete[] this->lpOutTies;
	delete[] this->lpInTies;

	this->lpOutTies = 0;
	this->lpInTies = 0;
}

/**
 * Destructs this network.
 */
Network::~Network() {
	fireNetworkDisposeEvent();
	this->deleteArrays();
}

// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the number of actors acting as tie senders in this network.
 */
int Network::n() const {
	return this->ln;
}

/**
 * Returns the number of actors acting as tie receivers in this network.
 */
int Network::m() const {
	return this->lm;
}

/**
 * Returns the total number of ties of this network
 */
int Network::tieCount() const {
	return this->ltieCount;
}

// ----------------------------------------------------------------------------
// Section: Basic structural operations
// ----------------------------------------------------------------------------

/**
 * This method sets the value of the tie from <i>i</i> to <i>j</i>.
 */
void Network::setTieValue(int i, int j, int v) {
	this->changeTieValue(i, j, v, REPLACE);
}

/**
 * This method increases the value of the tie from <i>i</i> to <i>j</i> by
 * <i>v</i> and returns the new value.
 */
int Network::increaseTieValue(int i, int j, int v) {
	return this->changeTieValue(i, j, v, INCREASE);
}

/**
 * This method changes the values of the tie from <i>i</i> to <i>j</i>
 * according to the specified type of change. The new value is returned as
 * the result.
 */
int Network::changeTieValue(int i, int j, int v, ChangeType type) {
	this->checkSenderRange(i);
	this->checkReceiverRange(j);

	// Retrieve the old value
	int oldValue = 0;
	std::map<int, int>& egoMap = lpOutTies[i];
	std::map<int, int>::iterator iter = egoMap.lower_bound(j);
	// we found the element
	if (iter != egoMap.end() && !egoMap.key_comp()(j, iter->first)) {
		oldValue = iter->second;
	}

	// Should we increase the value or replace?
	if (type == INCREASE) {
		v += oldValue;
	}

	// Update the maps of incoming and outgoing ties

	// if oldValue != 0 the (i,j) exists and we have to update the value
	// or remove the tie. Otherwise, we have to insert the tie (i,j) if
	// v!=0
	if (oldValue) {
		if (v == 0) {
			// A non-zero tie becomes 0. Just remove the corresponding
			// entries from the maps. Erasing an element pointed to by
			// an iterator is potentially faster than removing by key,
			// since we don't have to find the element.
			egoMap.erase(iter);
			this->lpInTies[j].erase(i);
		} else {
			// the value of the edge has been changed
			iter->second = v;
			this->lpInTies[j][i] = v;
		}
	} else if (v) {
		// iter points to the element right after j. Using this position
		// as a hint speeds things up.
		egoMap.insert(iter, std::map<int, int>::value_type(j, v));
		lpInTies[j].insert(std::map<int, int>::value_type(i, v));
	}
	// Remember that the network has changed
	this->lmodificationCount++;

	// Act on tie withdrawal or introduction.
	if (oldValue && !v) {
		// The (i,j) has been withdrawn.
		this->onTieWithdrawal(i, j);
	} else if (!oldValue && v) {
		// A new tie has been introduced.
		this->onTieIntroduction(i, j);
	}
	return v;
}

/**
 * Updates the state of this network to reflect the withdrawal of a tie
 * from actor <i>i</i> to actor <i>j</i>.
 */
void Network::onTieWithdrawal(int i, int j) {
	this->ltieCount--;
	// fire the withdrawal event
	fireWithdrawalEvent(i, j);
}

/**
 * Updates the state of this network to reflect the introduction of a tie
 * from actor <i>i</i> to actor <i>j</i>.
 */
void Network::onTieIntroduction(int i, int j) {
	this->ltieCount++;
	// fire the introduction even
	fireIntroductionEvent(i, j);
}

/**
 * Returns the value of the tie from <i>i</i> to <i>j</i>. The default is 0.
 */
int Network::tieValue(int i, int j) const {
	this->checkSenderRange(i);
	this->checkReceiverRange(j);

	// Look for the tie
	std::map<int, int>::const_iterator iter = this->lpOutTies[i].find(j);

	// The default value.
	int v = 0;

	if (iter != this->lpOutTies[i].end()) {
		// We have found an existing tie.
		v = iter->second;
	}

	return v;
}

/**
 * This method removes all ties from this network.
 */
void Network::clear() {
	// Clear the maps and reset the various degree counters.

	for (int i = 0; i < this->ln; i++) {
		this->lpOutTies[i].clear();
	}

	for (int i = 0; i < this->lm; i++) {
		this->lpInTies[i].clear();
	}

	// The ties are gone.
	this->ltieCount = 0;

	// The network has changed
	this->lmodificationCount++;

	// fire network clear event
	fireNetworkClearEvent();
}

/**
 * Removes all incoming ties of the given actor.
 */
void Network::clearInTies(int actor) {
	// We delegate to setTieValue such that various counters are updated
	// correctly.
	while (!this->lpInTies[actor].empty()) {
		int sender = this->lpInTies[actor].begin()->first;
		this->setTieValue(sender, actor, 0);
	}
}

/**
 * Removes all outgoing ties of the given actor.
 */
void Network::clearOutTies(int actor) {
	// We delegate to setTieValue such that various counters are updated
	// correctly.
	while (!this->lpOutTies[actor].empty()) {
		int receiver = this->lpOutTies[actor].begin()->first;
		this->setTieValue(actor, receiver, 0);
	}
}

// ----------------------------------------------------------------------------
// Section: Iterators
// ----------------------------------------------------------------------------

/**
 * Returns an iterator over all ties of this network.
 */
TieIterator Network::ties() const {
	return TieIterator(this);
}

/**
 * Returns an iterator over incoming ties of the actor <i>i</i>.
 */
IncidentTieIterator Network::inTies(int i) const {
	this->checkReceiverRange(i);
	return IncidentTieIterator(this->lpInTies[i]);
}

/**
 * Returns an iterator over outgoing ties of the actor <i>i</i> with
 * the receiver not less than the given bound.
 */
IncidentTieIterator Network::outTies(int i, int lowerBound) const {
	this->checkSenderRange(i);
	return IncidentTieIterator(this->lpOutTies[i], lowerBound);
}

/**
 * Returns an iterator over incoming ties of the actor <i>i</i> with
 * the sender not less than the given bound.
 */
IncidentTieIterator Network::inTies(int i, int lowerBound) const {
	this->checkReceiverRange(i);
	return IncidentTieIterator(this->lpInTies[i], lowerBound);
}

/**
 * Returns an iterator over outgoing ties of the actor <i>i</i>
 */
IncidentTieIterator Network::outTies(int i) const {
	this->checkSenderRange(i);
	return IncidentTieIterator(this->lpOutTies[i]);
}

// ----------------------------------------------------------------------------
// Section: Degrees
// ----------------------------------------------------------------------------

/**
 * Returns the number of incoming ties of the actor <i>i</i>.
 */
int Network::inDegree(int i) const {
	this->checkReceiverRange(i);
	return this->lpInTies[i].size();
}

/**
 * Returns the number of outgoing ties of the actor <i>i</i>.
 */
int Network::outDegree(int i) const {
	this->checkSenderRange(i);
	return this->lpOutTies[i].size();
}

// ----------------------------------------------------------------------------
// Section: Some useful statistics
// ----------------------------------------------------------------------------

/**
 * Returns the minimal tie value in this network. The implicit values
 * of 0 are also considered.
 */
int Network::minTieValue() const {
	// This method is linear in the total number of ties. Try maintaining
	// a sorted multi-set of tie values, if this method turns out to be
	// a bottleneck. It would add a log(ltieCount) overhead to the method
	// setTieValue(...), but make this method constant.

	int minValue = std::numeric_limits<int>::max();

	for (TieIterator iter = this->ties(); iter.valid(); iter.next()) {
		minValue = std::min(minValue, iter.value());
	}

	if (!this->complete()) {
		// The network is not complete, hence some tie variables are 0.
		minValue = std::min(minValue, 0);
	}

	return minValue;
}

/**
 * Returns the maximal tie value in this network. The implicit values
 * of 0 are also considered.
 */
int Network::maxTieValue() const {
	// This method is linear in the total number of ties. Try maintaining
	// a sorted multi-set of tie values, if this method turns out to be
	// a bottleneck. It would add a log(ltieCount) overhead to the method
	// setTieValue(...), but make this method constant.

	int maxValue = std::numeric_limits<int>::min();

	for (TieIterator iter = this->ties(); iter.valid(); iter.next()) {
		maxValue = std::max(maxValue, iter.value());
	}

	if (!this->complete()) {
		// The network is not complete, hence some tie variables are 0.
		maxValue = std::max(maxValue, 0);
	}

	return maxValue;
}

/**
 * This method returns the number of actors with ties to both <i>i</i> and
 * <i>j</i>.
 */
int Network::outTwoStarCount(int i, int j) const {
	this->checkReceiverRange(i);
	this->checkReceiverRange(j);
	return commonActorCount(this->inTies(i), this->inTies(j));
}

/**
 * This method returns the number of actors with ties from both <i>i</i> and
 * <i>j</i>.
 */
int Network::inTwoStarCount(int i, int j) const {
	this->checkSenderRange(i);
	this->checkSenderRange(j);
	return commonActorCount(this->outTies(i), this->outTies(j));
}

/**
 * Indicates that all ties are non-zero.
 */
bool Network::complete() const {
	// Return if all possible ties are non-zero.
	return this->ltieCount == this->maxTieCount();
}

/**
 * Tells whether there is an edge from <i>ego</i> to <i>alter</i>.
 */
bool Network::hasEdge(int ego, int alter) const {
	checkSenderRange(ego);
	checkReceiverRange(alter);
	return lpOutTies[ego].find(alter) != lpOutTies[ego].end();
}

/**
 * Returns the maximal possible number of ties in this network.
 */
int Network::maxTieCount() const {
	return this->ln * this->lm;
}

// ----------------------------------------------------------------------------
// Section: Helper checks
// ----------------------------------------------------------------------------

/**
 * Tests if the given actor is in the valid range of senders and throws an
 * std::out_of_range exception, if not.
 */
void Network::checkSenderRange(int i) const {
	if (i < 0 || i >= this->ln) {
		throw std::out_of_range(
				"The number " + toString(i) + " is not in the range [0,"
						+ toString(this->ln)
						+ ") of actors acting as senders of ties");
	}
}

/**
 * Tests if the given actor is in the valid range of receivers and throws an
 * std::out_of_range exception, if not.
 */
void Network::checkReceiverRange(int i) const {
	if (i < 0 || i >= this->lm) {
		throw std::out_of_range(
				"The number " + toString(i) + " is not in the range [0,"
						+ toString(this->lm)
						+ ") of actors acting as receivers of ties");
	}
}

/**
 * Adds the given <i>listener</i> from the network, if it is not yet attached.
 */
void Network::addNetworkChangeListener(
		INetworkChangeListener* const listener) {
	// ensure that the list is a set (no duplicates)
	std::list<INetworkChangeListener*>::iterator tmp = std::find(
			lNetworkChangeListener.begin(), lNetworkChangeListener.end(),
			listener);
	if (tmp == lNetworkChangeListener.end()) {
		lNetworkChangeListener.push_back(listener);
		listener->onInitializationEvent(*this);
	}
}

/**
 * Removes the given <i>listener</i> from the network.
 */
void Network::removeNetworkChangeListener(
		INetworkChangeListener* const listener) {
	std::list<INetworkChangeListener*>::iterator tmp = std::find(
			lNetworkChangeListener.begin(), lNetworkChangeListener.end(),
			listener);
	if (tmp != lNetworkChangeListener.end()) {
		// no need for erase remove since add ensures that the element is unique
		lNetworkChangeListener.erase(tmp);
	}
}

void Network::fireNetworkDisposeEvent() {
	for (std::list<INetworkChangeListener*>::const_iterator iter =
			lNetworkChangeListener.begin();
			iter != lNetworkChangeListener.end(); ++iter) {
		(*iter)->onNetworkDisposeEvent(*this);
	}
}

void Network::fireNetworkClearEvent() const {
	for (std::list<INetworkChangeListener*>::const_iterator iter =
			lNetworkChangeListener.begin();
			iter != lNetworkChangeListener.end(); ++iter) {
		(*iter)->onNetworkClearEvent(*this);
	}
}

/**
 * Inform all listeners that edge (ego,alter) has been inserted to the network.
 */
void Network::fireIntroductionEvent(int ego, int alter) const {
	for (std::list<INetworkChangeListener*>::const_iterator iter =
			lNetworkChangeListener.begin();
			iter != lNetworkChangeListener.end(); ++iter) {
		(*iter)->onTieIntroductionEvent(*this, ego, alter);
	}
}

/**
 * Tells whether the network is a one mode network or not.
 * @return <code>True</code> if the network is a one mode network,
 * <code>False</code> otherwise.
 */
bool Network::isOneMode() const {
	return false;
}

/**
 * Inform all listeners that edge (ego,alter) has been removed to the network.
 */
void Network::fireWithdrawalEvent(int ego, int alter) const {
	for (std::list<INetworkChangeListener*>::const_iterator iter =
			lNetworkChangeListener.begin();
			iter != lNetworkChangeListener.end(); ++iter) {
		(*iter)->onTieWithdrawalEvent(*this, ego, alter);
	}
}

}
