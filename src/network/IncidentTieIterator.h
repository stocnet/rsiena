/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IncidentTieIterator.h
 *
 * Description: This module defines the class IncidentTieIterator for
 * convenient iteration over incoming or outgoing ties of an actor.
 *****************************************************************************/

#ifndef INCIDENTTIEITERATOR_H_
#define INCIDENTTIEITERATOR_H_

#include <map>

#include "iterators/ITieIterator.h"

namespace siena {

/**
 * This class defines an iterator over incoming or outgoing ties of a specific
 * actor <i>i</i>. The ties are sorted in an increasing order of the neighbors
 * of <i>i</i>.
 */
class IncidentTieIterator: public ITieIterator {
	// The class Network needs access to the private constructor.
	friend class Network;
	friend class DistanceTwoLayer;

public:
	IncidentTieIterator();
	IncidentTieIterator(const IncidentTieIterator& rhs);


	int size() const {
		return -1;
	}
	/**
	 * Returns the neighbor incident to the current tie.
	 */
	inline int actor() const {
		if (valid()) {
			return lcurrent->first;
		}
		throw InvalidIteratorException();
	}

	/**
	 * Returns the value of the current tie.
	 */
	inline int value() const {
		if (valid()) {
			return lcurrent->second;
		}
		throw InvalidIteratorException();
	}

	/**
	 * Indicates if the iterator still points to a valid tie.
	 */
	inline bool valid() const {
		return lcurrent != lend;
	}

	/**
	 * Moves the iterator to the next tie.
	 */
	inline void next() {
		++lcurrent;
	}

	inline void reset() {
		lcurrent = lstart;
	}

	IncidentTieIterator* clone() const;

private:
	IncidentTieIterator(const std::map<int, int> & ties);
	IncidentTieIterator(const std::map<int, int> & ties, int lowerBound);


	/////////////////////////////////////////////////////////
	//NEVER CHANGE THIS ORDERING!!! CHECK INITIALIZATION LIST
	/////////////////////////////////////////////////////////

	// Points to the start element
	std::map<int, int>::const_iterator lstart;
	// Points to the current element in the underlying map
	std::map<int, int>::const_iterator lcurrent;

	// Points to the end of the underlying map
	std::map<int, int>::const_iterator lend;
};

}

#endif /*INCIDENTTIEITERATOR_H_*/
