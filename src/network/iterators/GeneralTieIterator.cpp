/*
 * GeneralTieIterator.cpp
 *
 *  Created on: May 13, 2015
 *      Author: ortmann
 */

#include <network/iterators/GeneralTieIterator.h>

namespace siena {

const Set_Operation Set_Operation::UNION = Set_Operation(0);
const Set_Operation Set_Operation::INTERSECTION = Set_Operation(1);
const Set_Operation Set_Operation::SYMDIFF = Set_Operation(2);
const Set_Operation Set_Operation::SET_MINUS = Set_Operation(3);
const Filter_Operation Filter_Operation::FILTER = Filter_Operation(0);
const Filter_Operation Filter_Operation::KEEP = Filter_Operation(1);

GeneralTieIterator::GeneralTieIterator() :
		ITieIterator(), //
		lPos(0), //
		rEntries(), //
		lSize(-1) {
}

GeneralTieIterator::GeneralTieIterator(ITieIterator& iter1, ITieIterator& iter2,
		Set_Operation opType) :
		ITieIterator(), //
		lPos(0) {
	init(iter1, iter2, opType);
	finalize(iter1, iter2);
}

GeneralTieIterator::GeneralTieIterator(ITieIterator& iter,
		const bool* const filter, const Filter_Operation opType) :
		ITieIterator(), //
		lPos(0) {
	if (opType == Filter_Operation::FILTER) {
		calcFilter(iter, filter, false);
	} else if (opType == Filter_Operation::KEEP) {
		calcFilter(iter, filter, true);
	} else {
		throw "unsupported operation type";
	}
	lSize = rEntries.size();
}

GeneralTieIterator::GeneralTieIterator(
		std::map<int, double>::const_iterator start,
		std::map<int, double>::const_iterator end) :
		ITieIterator(), //
		lPos(0) {
	while (start != end) {
		rEntries.push_back(start->first);
		++start;
	}
	lSize = rEntries.size();
}

GeneralTieIterator::GeneralTieIterator(const int actor) :
		ITieIterator(), //
		lPos(0), //
		lSize(1) {
	rEntries.push_back(actor);
}

GeneralTieIterator::GeneralTieIterator(std::vector<int>::const_iterator start,
		std::vector<int>::const_iterator end) :
		ITieIterator(), //
		lPos(0), //
		rEntries(start, end), //
		lSize(rEntries.size()) {
}

GeneralTieIterator::~GeneralTieIterator() {
}

GeneralTieIterator* GeneralTieIterator::clone() const {
	return new GeneralTieIterator(*this);
}

GeneralTieIterator::GeneralTieIterator(const GeneralTieIterator& rhs) :
		ITieIterator(), //
		lPos(rhs.lPos), //
		rEntries(rhs.rEntries), //
		lSize(rhs.lSize) {
}

void GeneralTieIterator::init(ITieIterator& iter1, ITieIterator& iter2,
		const Set_Operation opType) {
	if (opType == Set_Operation::UNION) {
		calcUnion(iter1, iter2);
	} else if (opType == Set_Operation::INTERSECTION) {
		calcIntersection(iter1, iter2);
	} else if (opType == Set_Operation::SET_MINUS) {
		calcSetMinus(iter1, iter2);
	} else {
		throw "no such set operation implemented";
	}
}

void GeneralTieIterator::finalize(ITieIterator& iter1, ITieIterator& iter2) {
	lSize = rEntries.size();
	iter1.reset();
	iter2.reset();
}

void GeneralTieIterator::calcUnion(ITieIterator& iter1, ITieIterator& iter2) {
	int actorIter1;
	int actorIter2;
// calc the union
	while (iter1.valid() && iter2.valid()) {
		actorIter1 = iter1.actor();
		actorIter2 = iter2.actor();
		if (actorIter1 < actorIter2) {
			rEntries.push_back(actorIter1);
			iter1.next();
		} else if (actorIter1 > actorIter2) {
			rEntries.push_back(actorIter2);
			iter2.next();
		} else {
			rEntries.push_back(actorIter1);
			iter1.next();
			iter2.next();
		}
	}
// add remaining entries
	if (iter1.valid()) {
		while (iter1.valid()) {
			rEntries.push_back(iter1.actor());
			iter1.next();
		}
	}
	if (iter2.valid()) {
		while (iter2.valid()) {
			rEntries.push_back(iter2.actor());
			iter2.next();
		}
	}
}

void GeneralTieIterator::calcIntersection(ITieIterator& iter1,
		ITieIterator& iter2) {
	int actorIter1;
	int actorIter2;
// calc the union
	while (iter1.valid() && iter2.valid()) {
		actorIter1 = iter1.actor();
		actorIter2 = iter2.actor();
		if (actorIter1 < actorIter2) {
			iter1.next();
		} else if (actorIter1 > actorIter2) {
			iter2.next();
		} else {
			rEntries.push_back(actorIter1);
			iter1.next();
			iter2.next();
		}
	}
}

void GeneralTieIterator::calcSymDiff(ITieIterator& iter1, ITieIterator& iter2) {
	int actorIter1;
	int actorIter2;
// calc the union
	while (iter1.valid() && iter2.valid()) {
		actorIter1 = iter1.actor();
		actorIter2 = iter2.actor();
		if (actorIter1 < actorIter2) {
			rEntries.push_back(actorIter1);
			iter1.next();
		} else if (actorIter1 > actorIter2) {
			rEntries.push_back(actorIter2);
			iter2.next();
		} else {
			iter1.next();
			iter2.next();
		}
	}
// add remaining entries
	if (iter1.valid()) {
		while (iter1.valid()) {
			rEntries.push_back(iter1.actor());
			iter1.next();
		}
	}
	if (iter2.valid()) {
		while (iter2.valid()) {
			rEntries.push_back(iter2.actor());
			iter2.next();
		}
	}
}

void GeneralTieIterator::calcSetMinus(ITieIterator& iter1,
		ITieIterator& iter2) {
	int actorIter1;
	int actorIter2;
// calc the union
	while (iter1.valid() && iter2.valid()) {
		actorIter1 = iter1.actor();
		actorIter2 = iter2.actor();
		if (actorIter1 < actorIter2) {
			while (iter1.valid() && iter1.actor() < actorIter2) {
				rEntries.push_back(iter1.actor());
				iter1.next();
			}
			if (!iter1.valid()) {
				return;
			}
		}
		actorIter1 = iter1.actor();
		if (actorIter2 < actorIter1) {
			while (iter2.valid() && iter2.actor() < actorIter1) {
				iter2.next();
			}
			if (!iter2.valid()) {
				break;
			}
		}
		actorIter2 = iter2.actor();
		if (actorIter1 == actorIter2) {
			iter1.next();
			iter2.next();
		}
	}
	while (iter1.valid()) {
		rEntries.push_back(iter1.actor());
		iter1.next();
	}
}

void GeneralTieIterator::calcFilter(ITieIterator& iter,
		const bool* const filter, const bool keep) {
	while (iter.valid()) {
		if (filter[iter.actor()] == keep) {
			rEntries.push_back(iter.actor());
		}
		iter.next();
	}
	iter.reset();
}

} /* namespace siena */
