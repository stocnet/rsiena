#include "CommonNeighborIterator.h"
#include "../utils/Utils.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: Construction area
// ----------------------------------------------------------------------------

/**
 * Creates an iterator over actors common to both of the given incident
 * tie iterators.
 */
CommonNeighborIterator::CommonNeighborIterator(const IncidentTieIterator& iter1,
		const IncidentTieIterator& iter2) :
		ITieIterator(), //
		liter1(iter1), //
		liter2(iter2) {
	// Make sure they point to the first common actor.
	this->skipMismatches();
}

// ----------------------------------------------------------------------------
// Section: Public interface
// ----------------------------------------------------------------------------

/**
 * Indicates if there are still some common actors to be reported.
 */
bool CommonNeighborIterator::valid() const {
	return this->liter1.valid() && this->liter2.valid();
}

/**
 * Returns the current common actor.
 */
int CommonNeighborIterator::actor() const {
	if (!this->valid()) {
		throw InvalidIteratorException();
	}

	// Both iterators point to the same actor, so we can use any iterator.
	return this->liter1.actor();
}

/**
 * Moves on to the next common actor if any.
 */
void CommonNeighborIterator::next() {
	if (!this->valid()) {
		throw InvalidIteratorException();
	}

	// Advance both iterators until they point to the next common actor
	// or we run out of actors.

	this->liter1.next();
	this->liter2.next();
	this->skipMismatches();
}

CommonNeighborIterator* CommonNeighborIterator::clone() const {
	return new CommonNeighborIterator(*this);
}

// ----------------------------------------------------------------------------
// Section: Private methods
// ----------------------------------------------------------------------------

/**
 * Advances both iterators until they point to the same actor. If there is
 * no such an actor, one of the iterators becomes invalid and we stop.
 */
void CommonNeighborIterator::skipMismatches() {
	while (this->liter1.valid() && this->liter2.valid()
			&& this->liter1.actor() != this->liter2.actor()) {
		while (liter1.valid() && liter1.actor() < liter2.actor()) {
			liter1.next();
		}
		if (!liter1.valid()) {
			return;
		}
		while (liter2.valid() && liter2.actor() < liter1.actor()) {
			liter2.next();
		}
	}
}

CommonNeighborIterator::CommonNeighborIterator(
		const CommonNeighborIterator& rhs) :
		ITieIterator(rhs), //
		liter1(rhs.liter1), //
		liter2(rhs.liter2) {
	// note there is no need to skip matches this has been down by rhs
}

}

void siena::CommonNeighborIterator::reset() {
	liter1.reset();
	liter2.reset();
}
