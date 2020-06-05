/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: UnionTieIterator.h
 *
 * Description: This module defines an iterator that only the union
 * of two iterators.
 *****************************************************************************/

#ifndef UNIONTIEITERATOR_H_
#define UNIONTIEITERATOR_H_

#include "GeneralTieIterator.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: UnionTieIterator class
// ----------------------------------------------------------------------------

class UnionTieIterator: public GeneralTieIterator {
public:

	UnionTieIterator(ITieIterator& iter1, ITieIterator& iter2) :
			GeneralTieIterator(iter1, iter2, Set_Operation::UNION) {
	}

	virtual ~UnionTieIterator() {
	}

	UnionTieIterator * clone() const {
		return new UnionTieIterator(*this);
	}

protected:
	UnionTieIterator(const UnionTieIterator& rhs) :
			GeneralTieIterator(rhs) {
	}

};

} /* namespace siena */

#endif /* UNIONTIEITERATOR_H_ */
