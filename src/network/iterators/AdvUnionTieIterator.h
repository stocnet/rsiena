/*
 * AdvUnionTieIterator.h
 *
 *  Created on: May 13, 2015
 *      Author: ortmann
 */

#ifndef ADVUNIONTIEITERATOR_H_
#define ADVUNIONTIEITERATOR_H_

#include "GeneralTieIterator.h"

namespace siena {

class AdvUnionTieIterator: public GeneralTieIterator {
public:
	AdvUnionTieIterator(int iter1ID, int iter2ID, ITieIterator& iter1,
			ITieIterator& iter2);
	virtual ~AdvUnionTieIterator();

	bool isCommon() const;

	int getInactiveIterID() const;

	virtual AdvUnionTieIterator * clone() const;

protected:

	AdvUnionTieIterator(const AdvUnionTieIterator& rhs);

	void common(const int actor, const bool isCommon, const int iterID);

private:

	int literOneID{};

	int literTwoID{};

	std::vector<bool> rCommon;

	std::vector<int> rInactiveIter;

	void calcAdvUnion(ITieIterator& iter1,ITieIterator& iter2);
};

} /* namespace siena */

#endif /* ADVUNIONTIEITERATOR_H_ */
