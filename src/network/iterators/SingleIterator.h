/*
 * SingleIterator.h
 *
 *  Created on: 24.06.2014
 *      Author: ortmann
 */

#ifndef SINGLEITERATOR_H_
#define SINGLEITERATOR_H_

#include "GeneralTieIterator.h"

namespace siena {

class SingleIterator: public GeneralTieIterator {
public:
	SingleIterator(int actor) :
			GeneralTieIterator(actor) {
	}

	virtual ~SingleIterator(){
	}

	virtual SingleIterator * clone() const {
		return new SingleIterator(*this);
	}

protected:

	SingleIterator(const SingleIterator& rhs) :
			GeneralTieIterator(rhs) {
	}
};

} /* namespace siena */
#endif /* SINGLEITERATOR_H_ */
