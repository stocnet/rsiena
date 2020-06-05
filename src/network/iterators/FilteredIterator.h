/*
 * FilteredIterator.h
 *
 *  Created on: 24.06.2014
 *      Author: ortmann
 */

#ifndef FILTEREDITERATOR_H_
#define FILTEREDITERATOR_H_

#include "GeneralTieIterator.h"

namespace siena {

class FilteredIterator: public GeneralTieIterator {
public:
	FilteredIterator(ITieIterator& iter1, const bool* const filter) :
			GeneralTieIterator(iter1, filter, Filter_Operation::KEEP) {
	}

	virtual ~FilteredIterator(){

	}

	FilteredIterator * clone() const {
		return new FilteredIterator(*this);
	}

protected:

	FilteredIterator(const FilteredIterator& rhs) :
			GeneralTieIterator(rhs) {
	}

};

} /* namespace siena */
#endif /* FILTEREDITERATOR_H_ */
