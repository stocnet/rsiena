/*
 * AdvUnionTieIterator.cpp
 *
 *  Created on: May 13, 2015
 *      Author: ortmann
 */

#include "network/iterators/AdvUnionTieIterator.h"

namespace siena {

AdvUnionTieIterator::AdvUnionTieIterator(int iter1ID, int iter2ID,
		ITieIterator& iter1, ITieIterator& iter2) :
		GeneralTieIterator(), //
		literOneID(iter1ID), //
		literTwoID(iter2ID) {
	calcAdvUnion(iter1, iter2);
	finalize(iter1, iter2);
}

AdvUnionTieIterator::~AdvUnionTieIterator() {
}

bool AdvUnionTieIterator::isCommon() const {
	return rCommon[lPos];
}

int AdvUnionTieIterator::getInactiveIterID() const {
	return rInactiveIter[lPos];
}

AdvUnionTieIterator* AdvUnionTieIterator::clone() const {
	return new AdvUnionTieIterator(*this);
}

AdvUnionTieIterator::AdvUnionTieIterator(const AdvUnionTieIterator& rhs) :
		GeneralTieIterator(rhs), //
		literOneID(rhs.literOneID), //
		literTwoID(rhs.literTwoID), //
		rCommon(rhs.rCommon), //
		rInactiveIter(rhs.rInactiveIter) {
}

void AdvUnionTieIterator::common(const int actor, const bool isCommon,
		const int iterID) {
	rCommon.push_back(isCommon);
	if (iterID == 1) {
		rInactiveIter.push_back(literTwoID);
	} else if (iterID == 2) {
		rInactiveIter.push_back(literOneID);
	} else {
		rInactiveIter.push_back(-1);
	}
}

void AdvUnionTieIterator::calcAdvUnion(ITieIterator& iter1,
		ITieIterator& iter2) {
	int actorIter1;
	int actorIter2;
// calc the union
	while (iter1.valid() && iter2.valid()) {
		actorIter1 = iter1.actor();
		actorIter2 = iter2.actor();
		if (actorIter1 < actorIter2) {
			rCommon.push_back(false);
			rInactiveIter.push_back(literTwoID);
			rEntries.push_back(actorIter1);
			iter1.next();
		} else if (actorIter1 > actorIter2) {
			rCommon.push_back(false);
			rInactiveIter.push_back(literOneID);
			rEntries.push_back(actorIter2);
			iter2.next();
		} else {
			rCommon.push_back(true);
			rInactiveIter.push_back(literOneID);
			rEntries.push_back(actorIter1);
			iter1.next();
			iter2.next();
		}
	}
// add remaining entries
	if (iter1.valid()) {
		while (iter1.valid()) {
			rCommon.push_back(false);
			rInactiveIter.push_back(literTwoID);
			rEntries.push_back(iter1.actor());
			iter1.next();
		}
	}
	if (iter2.valid()) {
		while (iter2.valid()) {
			rCommon.push_back(false);
			rInactiveIter.push_back(literOneID);
			rEntries.push_back(iter2.actor());
			iter2.next();
		}
	}
}

} /* namespace siena */
