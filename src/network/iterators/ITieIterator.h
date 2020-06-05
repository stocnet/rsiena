/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ITieIterator.h
 *
 * Description: This interface defines the common operations supported by
 * any TieIterator.
 *****************************************************************************/

#ifndef ITIEITERATOR_H_
#define ITIEITERATOR_H_

#include "../../utils/Utils.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: ITieIterator interface
// ----------------------------------------------------------------------------

class ITieIterator {
public:

	/**
	 * Destructor.
	 */
	virtual ~ITieIterator() {
	}

	/**
	 * Moves the iterator to the next position.
	 */
	virtual void next() = 0;

	/**
	 * Returns the actor at the current position. Note that the
	 * behavior of this function is undefined if the iterator is
	 * invalid.
	 * @return The actor at the current position.
	 */
	virtual int actor() const = 0;

	/**
	 * Tells whether the current position is valid or not.
	 * @return <code>True</code> indicating that the current position
	 * is valid and <code>False</code> otherwise.
	 */
	virtual bool valid() const = 0;


	/**
	 * Sets the iterator to the beginning
	 */
	virtual void reset() = 0;

	/**
	 * Returns the size of the iterator
	 */
	virtual int size() const = 0 ;


	/**
	 * Creates an identical copy of the iterator.
	 */
	virtual ITieIterator* clone() const = 0;
protected:

	/**
	 * Constructor.
	 */
	ITieIterator() {
	}

	/**
	 * Copy assignment Constructor.
	 */
	ITieIterator(const ITieIterator&) {
	}

	/**
	 * Assignment operator.
	 */
	ITieIterator& operator=(const ITieIterator&) {
		return *this;
	}
};

} /* namespace siena */
#endif /* ITIEITERATOR_H_ */
