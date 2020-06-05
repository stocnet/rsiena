/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OneModeNetwork.cpp
 *
 * Description: This file contains the implementation of the
 * OneModeNetwork class.
 *****************************************************************************/

#include <limits>
#include <stdexcept>
#include <string>

#include "OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "network/INetworkChangeListener.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: Construction and destruction
// ----------------------------------------------------------------------------

/**
 * Constructs a one-mode network.
 * @param[in] n the number of actors in the network
 * @param[in] loopsPermitted indicates if ties with equal senders and receivers
 * are permitted
 */
OneModeNetwork::OneModeNetwork(int n, bool loopsPermitted) :
		Network(n, n) {
	this->lloopsPermitted = loopsPermitted;

	// Initialize the reciprocal degree counters

	this->lpReciprocalDegree = new int[n];

	for (int i = 0; i < n; i++) {
		this->lpReciprocalDegree[i] = 0;
	}
}

/**
 * Constructs a copy of the given one-mode network.
 */
OneModeNetwork::OneModeNetwork(const OneModeNetwork & rNetwork) :
		Network(rNetwork) {
	// Copy the fields
	this->lloopsPermitted = rNetwork.lloopsPermitted;

	// Copy the reciprocal degree counters

	this->lpReciprocalDegree = new int[rNetwork.n()];

	for (int i = 0; i < rNetwork.n(); i++) {
		this->lpReciprocalDegree[i] = rNetwork.lpReciprocalDegree[i];
	}
}

/**
 * Assigns the contents of the given network to this network.
 */
OneModeNetwork & OneModeNetwork::operator=(const OneModeNetwork & rNetwork) {
	if (this != &rNetwork) {
		// Let the base class do its part
		Network::operator=(rNetwork);

		// Now copy our own fields
		this->lloopsPermitted = rNetwork.lloopsPermitted;

		// Reallocate the reciprocal degree counters

		delete[] this->lpReciprocalDegree;
		this->lpReciprocalDegree = new int[this->n()];

		// Copy the reciprocal degree counters

		for (int i = 0; i < rNetwork.n(); i++) {
			this->lpReciprocalDegree[i] = rNetwork.lpReciprocalDegree[i];
		}
	}
	for (std::list<INetworkChangeListener*>::const_iterator iter =
			lNetworkChangeListener.begin();
			iter != lNetworkChangeListener.end(); ++iter) {
		(*iter)->onInitializationEvent(*this);
	}

	return *this;
}

Network * OneModeNetwork::clone() const {
	return new OneModeNetwork(*this);
}

/**
 * Deallocates this network.
 */
OneModeNetwork::~OneModeNetwork() {
	delete[] this->lpReciprocalDegree;
	this->lpReciprocalDegree = 0;
}

// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Indicates if loops are permitted.
 */
bool OneModeNetwork::loopsPermitted() const {
	return this->lloopsPermitted;
}

// ----------------------------------------------------------------------------
// Section: Basic structural operations
// ----------------------------------------------------------------------------

/**
 * This method changes the values of the tie from <i>i</i> to <i>j</i>
 * according to the specified type of change. The new value is returned as
 * the result.
 */
int OneModeNetwork::changeTieValue(int i, int j, int v, ChangeType type) {
	if (i == j && !this->lloopsPermitted) {
		throw std::invalid_argument("Loops are not permitted for this network");
	}

	return Network::changeTieValue(i, j, v, type);
}

/**
 * Updates the state of this network to reflect the withdrawal of a tie
 * from actor <i>i</i> to actor <i>j</i>.
 */
void OneModeNetwork::onTieWithdrawal(int i, int j) {
	// Call the base method first
	Network::onTieWithdrawal(i, j);

	// A reciprocal tie might be lost.

	if (i == j) {
		// Careful with loops!
		this->lpReciprocalDegree[i]--;
	} else if (this->tieValue(j, i)) {
		this->lpReciprocalDegree[i]--;
		this->lpReciprocalDegree[j]--;
	}
}

/**
 * Updates the state of this network to reflect the introduction of a tie
 * from actor <i>i</i> to actor <i>j</i>.
 */
void OneModeNetwork::onTieIntroduction(int i, int j) {
	Network::onTieIntroduction(i, j);

	// The tie might be reciprocated.

	if (i == j) {
		// Careful with loops!
		this->lpReciprocalDegree[i]++;
	} else if (this->tieValue(j, i)) {
		this->lpReciprocalDegree[i]++;
		this->lpReciprocalDegree[j]++;
	}
}

/**
 * This method removes all ties from this network.
 */
void OneModeNetwork::clear() {
	// Let the base class do its part
	Network::clear();

	// Reset the degree counters.

	for (int i = 0; i < this->n(); i++) {
		this->lpReciprocalDegree[i] = 0;
	}
}

// ----------------------------------------------------------------------------
// Section: Iterators
// ----------------------------------------------------------------------------

/**
 * Returns an iterator over reciprocated ties of the actor <i>i</i>.
 */
CommonNeighborIterator OneModeNetwork::reciprocatedTies(int i) const {
	this->checkSenderRange(i);
	return CommonNeighborIterator(this->inTies(i), this->outTies(i));
}

/**
 * Returns an iterator over reciprocated ties of the actor <i>i</i> with
 * the alter not less than the given bound.
 */
CommonNeighborIterator OneModeNetwork::reciprocatedTies(int i,
		int lowerBound) const {
	this->checkSenderRange(i);
	return CommonNeighborIterator(this->inTies(i, lowerBound),
			this->outTies(i, lowerBound));
}

// ----------------------------------------------------------------------------
// Section: Degrees
// ----------------------------------------------------------------------------

/**
 * Returns the number of reciprocated ties of the actor <i>i</i>.
 */
int OneModeNetwork::reciprocalDegree(int i) const {
	this->checkSenderRange(i);
	return this->lpReciprocalDegree[i];
}

// ----------------------------------------------------------------------------
// Section: Some useful statistics and properties
// ----------------------------------------------------------------------------

/**
 * Indicates if all ties are reciprocated with the same value.
 */
bool OneModeNetwork::symmetric() const {
	// The current implementation is linear in the total number of ties.
	// The time complexity can be reduced to a constant by maintaining
	// the number of non-symmetric ties.

	// Assume the network is symmetric until we can disprove it.
	bool rc = true;

	// Test the incoming and outgoing ties of each actor in turn.

	for (int i = 0; i < this->n() && rc; i++) {
		if (this->outDegree(i) == this->inDegree(i)) {
			IncidentTieIterator outIter = this->outTies(i);
			IncidentTieIterator inIter = this->inTies(i);

			// No need to test both iterators for validity, as the numbers
			// of incoming and outgoing ties are the same.

			while (outIter.valid() && rc) {
				if (outIter.actor() != inIter.actor()
						|| outIter.value() != inIter.value()) {
					// Found a mismatch.
					rc = false;
				}

				outIter.next();
				inIter.next();
			}
		} else {
			// The numbers of incoming and outgoing ties differ, which
			// destroys the symmetry immediately.

			rc = false;
		}
	}

	return rc;
}

/**
 * Returns the number of two-paths from <i>i</i> to <i>j</i>.
 */
int OneModeNetwork::twoPathCount(int i, int j) const {
	return this->truncatedTwoPathCount(i, j, std::numeric_limits<int>::max());
}

/**
 * Returns the number of two-paths from <i>i</i> to <i>j</i> truncated at the
 * given threshold value.
 */
int OneModeNetwork::truncatedTwoPathCount(int i, int j, int threshold) const {
	this->checkSenderRange(i);
	this->checkReceiverRange(j);

	// Iterate the outgoing ties of i and incoming ties of j simultaneously
	// and count the number of matching neighbors. Stop as soon as the
	// threshold value is reached.

	IncidentTieIterator outIter = this->outTies(i);
	IncidentTieIterator inIter = this->inTies(j);
	int count = 0;

	while (outIter.valid() && inIter.valid() && count < threshold) {
		if (outIter.actor() < inIter.actor()) {
			outIter.next();
		} else if (outIter.actor() > inIter.actor()) {
			inIter.next();
		} else {
			count++;
			outIter.next();
			inIter.next();
		}
	}

	return count;
}

/**
 * This method indicated that there are no two-paths from <i>i</i> to <i>j</i>
 * via an intermediate actor below the given bound.
 */
bool OneModeNetwork::noTwoPaths(int i, int j,
		int intermediateActorUpperBound) const {
	this->checkSenderRange(i);
	this->checkReceiverRange(j);

	// Iterate the outgoing ties of i and incoming ties of j simultaneously
	// and check if there are no matching neighbors below the given bound.
	// We can terminate the loop as soon as one of the iterators reaches
	// the bound.

	IncidentTieIterator outIter = this->outTies(i);
	IncidentTieIterator inIter = this->inTies(j);
	bool found = false;

	while (outIter.valid() && inIter.valid() && !found
			&& outIter.actor() < intermediateActorUpperBound
			&& inIter.actor() < intermediateActorUpperBound) {
		if (outIter.actor() < inIter.actor()) {
			outIter.next();
		} else if (outIter.actor() > inIter.actor()) {
			inIter.next();
		} else {
			found = true;
		}
	}

	return !found;
}

/**
 * This method indicated that there is at least one two-path from <i>i</i>
 * to <i>j</i>.
 */
bool OneModeNetwork::existsTwoPath(int i, int j) const {
	return !this->noTwoPaths(i, j, std::numeric_limits<int>::max());
}

/**
 * This method tests if there are at most <i>k</i> two-paths from <i>i</i>
 * to <i>j</i>.
 * If the test is positive, the number of such paths is stored
 * in the variable <i>twoPathCount</i>.
 */
bool OneModeNetwork::atMostKTwoPaths(int i, int j, int k,
		int & twoPathCount) const {
	this->checkSenderRange(i);
	this->checkReceiverRange(j);

	// Iterate the outgoing ties of i and incoming ties of j simultaneously
	// and count the number of matching neighbors. Stop as soon as their number
	// exceeds k.

	IncidentTieIterator outIter = this->outTies(i);
	IncidentTieIterator inIter = this->inTies(j);
	twoPathCount = 0;

	while (outIter.valid() && inIter.valid() && twoPathCount <= k) {
		if (outIter.actor() < inIter.actor()) {
			outIter.next();
		} else if (outIter.actor() > inIter.actor()) {
			inIter.next();
		} else {
			twoPathCount++;
			outIter.next();
			inIter.next();
		}
	}

	return twoPathCount <= k;
}

/**
 * This method counts the number of actors participating in three or all four
 * ties with actors <i>i</i> and <i>j</i>. These numbers are stored in
 * <i>n3</i> and <i>n4</i>, respectively.
 */
void OneModeNetwork::neighborCensus(int i, int j, int & n3, int & n4) const {
	this->checkSenderRange(i);
	this->checkSenderRange(j);

	// There are four relevant iterators for actors i and j -- two for
	// incoming ties and two for outgoing ties.

	IncidentTieIterator iterators[] = { this->inTies(i), this->outTies(i),
			this->inTies(j), this->outTies(j) };

	// How many of these iterators are valid?

	int validIteratorCount = 0;

	for (int i = 0; i < 4; i++) {
		if (iterators[i].valid()) {
			validIteratorCount++;
		}
	}

	// Initialize the counters.

	n3 = 0;
	n4 = 0;

	// While there are at least three valid iterators, there is a chance
	// to encounter an actor involved in at least three ties with i and j.

	while (validIteratorCount >= 3) {
		// Find the minimum of the current actors of the valid iterators.

		int minActor = std::numeric_limits<int>::max();

		for (int i = 0; i < 4; i++) {
			if (iterators[i].valid()) {
				minActor = std::min(minActor, iterators[i].actor());
			}
		}

		// How many iterators have this minimum as the current actor?
		int count = 0;

		for (int i = 0; i < 4; i++) {
			if (iterators[i].valid() && iterators[i].actor() == minActor) {
				count++;

				// Advance the iterator right away and
				// check if it is still valid

				iterators[i].next();

				if (!iterators[i].valid()) {
					validIteratorCount--;
				}
			}
		}

		// Update the counters.

		if (count == 3) {
			n3++;
		} else if (count == 4) {
			n4++;
		}
	}
}

/**
 * @copydoc Network::isOneMode()
 */
bool OneModeNetwork::isOneMode() const {
	return true;
}

/**
 * Returns the maximal possible number of ties in this network.
 */
int OneModeNetwork::maxTieCount() const {
	int count = this->n() * this->n();

	if (!this->lloopsPermitted) {
		count = this->n() * (this->n() - 1);
	}

	return count;
}

}
