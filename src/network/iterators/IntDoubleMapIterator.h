/*
 * IntDoubleMapIterator.h
 *
 *  Created on: 25.06.2014
 *      Author: ortmann
 */

#ifndef INTDOUBLEMAPITERATOR_H_
#define INTDOUBLEMAPITERATOR_H_

#include "GeneralTieIterator.h"
#include <map>

namespace siena {

class IntDoubleMapIterator: public GeneralTieIterator {
public:
	IntDoubleMapIterator(std::map<int, double>::const_iterator start,
			std::map<int, double>::const_iterator end) :
			GeneralTieIterator(start, end) {
	}

	virtual ~IntDoubleMapIterator(){

	}

	virtual IntDoubleMapIterator * clone() const {
		return new IntDoubleMapIterator(*this);
	}

protected:

	IntDoubleMapIterator(const IntDoubleMapIterator& rhs) :
			GeneralTieIterator(rhs) {
	}
};

} /* namespace siena */
#endif /* INTDOUBLEMAPITERATOR_H_ */
