#include "UnionNeighborIterator.h"
#include "../utils/Utils.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: Construction area
// ----------------------------------------------------------------------------

/**
 * Creates an iterator over actors that are elements of either of the given incident
 * tie iterators.
 */
UnionNeighborIterator::UnionNeighborIterator(const IncidentTieIterator& iter1,
		const IncidentTieIterator& iter2) :
		ITieIterator(), //
		liter1(iter1), //
		liter2(iter2) {

}

// ----------------------------------------------------------------------------
// Section: Public interface
// ----------------------------------------------------------------------------

/**
 * Indicates if there are still some common actors to be reported.
 */
bool UnionNeighborIterator::valid() const {
	return (this->liter1.valid() || this->liter2.valid()); // changed to OR
}

/**
 * Returns the current common actor.
 */
int UnionNeighborIterator::actor() const {
	if (!this->valid()) {
		throw InvalidIteratorException();
	}

	// return the minimum actor of the valid pointers
	if (this->liter1.valid() && this->liter2.valid())
		return std::min(this->liter1.actor(), this->liter2.actor());
	if (this->liter1.valid())
		return(this->liter1.actor());
	return this->liter2.actor();
}

/**
 * Moves on to the next common actor if any.
 */
void UnionNeighborIterator::next() {
	if (!this->valid()) {
		throw InvalidIteratorException();
	}

	// increase the iterator that points at the minium value
	// if both point at the same value, increase both
	if (this->liter1.valid() && this->liter2.valid()){
		if (this->liter1.actor() == this->liter2.actor()){
			this->liter1.next();
			this->liter2.next();
			return;
		}
		if (this->liter1.actor() < this->liter2.actor())
			this->liter1.next();
		else this->liter2.next();
		return;
	}

	// if only one iterator is valid, increase that one
	if (this->liter1.valid()) this->liter1.next();
	if (this->liter2.valid()) this->liter2.next();

}

UnionNeighborIterator* UnionNeighborIterator::clone() const {
	return new UnionNeighborIterator(*this);
}

// ----------------------------------------------------------------------------
// Section: Private methods
// ----------------------------------------------------------------------------


UnionNeighborIterator::UnionNeighborIterator(
		const UnionNeighborIterator& rhs) :
		ITieIterator(rhs), //
		liter1(rhs.liter1), //
		liter2(rhs.liter2) {
	// note there is no need to skip matches this has been down by rhs
}

void siena::UnionNeighborIterator::reset() {
  liter1.reset();
  liter2.reset();
}


}
