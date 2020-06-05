/*
 * IntVecIterator.h
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#ifndef INTVECITERATOR_H_
#define INTVECITERATOR_H_

#include "GeneralTieIterator.h"

namespace siena {

class IntVecIterator: public GeneralTieIterator {
public:
	IntVecIterator(std::vector<int>::const_iterator start,
			std::vector<int>::const_iterator end) :
			GeneralTieIterator(start, end) {
	}

	virtual ~IntVecIterator(){
	}

	virtual IntVecIterator * clone() const {
		return new IntVecIterator(*this);
	}

protected:

	IntVecIterator(const IntVecIterator& rhs) :
			GeneralTieIterator(rhs) {
	}

};

} /* namespace siena */
#endif /* INTVECITERATOR_H_ */
