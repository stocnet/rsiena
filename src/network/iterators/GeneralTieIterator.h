/*
 * GeneralTieIterator.h
 *
 *  Created on: May 13, 2015
 *      Author: ortmann
 */

#ifndef GENERALTIEITERATOR_H_
#define GENERALTIEITERATOR_H_

#include "ITieIterator.h"

#include <map>
#include <vector>

namespace siena {

class Set_Operation {
public:
	static const Set_Operation UNION;
	static const Set_Operation INTERSECTION;
	static const Set_Operation SYMDIFF;
	static const Set_Operation SET_MINUS;

	bool operator==(const Set_Operation& rhs) const {
		return value == rhs.value;
	}

private:
	Set_Operation(int val) :
			value(val) {
	}
	int value {};
};

class Filter_Operation {
public:
	static const Filter_Operation FILTER;
	static const Filter_Operation KEEP;

	bool operator==(const Filter_Operation& rhs) const {
		return value == rhs.value;
	}

private:
	Filter_Operation(int val) :
			value(val) {
	}
	int value {};
};

class GeneralTieIterator: public ITieIterator {
public:
	GeneralTieIterator(ITieIterator& iter1, ITieIterator& iter2,
			const Set_Operation opType);
	GeneralTieIterator(ITieIterator& iter, const bool* const filter,
			const Filter_Operation opType);
	GeneralTieIterator(std::vector<int>::const_iterator start,
			std::vector<int>::const_iterator end);
	GeneralTieIterator(std::map<int, double>::const_iterator start,
			std::map<int, double>::const_iterator end);
	GeneralTieIterator(const int actor);
	virtual ~GeneralTieIterator();

	inline void next() {
		++lPos;
	}

	inline int actor() const {
		if (valid()) return rEntries[lPos];
		return -1;
	}

	inline bool valid() const {
		return lPos < lSize;
	}

	inline void reset() {
		lPos = 0;
	}

	inline int size() const {
		return lSize;
	}

	virtual GeneralTieIterator * clone() const;

protected:
	GeneralTieIterator();
	GeneralTieIterator(const GeneralTieIterator& rhs);
	int lPos;
	std::vector<int> rEntries;
	int lSize;
	void finalize(ITieIterator& iter1, ITieIterator& iter2);

private:
	void init(ITieIterator& iter1, ITieIterator& iter2,
			const Set_Operation opType);
	void calcUnion(ITieIterator& iter1, ITieIterator& iter2);
	void calcIntersection(ITieIterator& iter1, ITieIterator& iter2);
	void calcSymDiff(ITieIterator& iter1, ITieIterator& iter2);
	void calcSetMinus(ITieIterator& iter1, ITieIterator& iter2);
	void calcFilter(ITieIterator& iter, const bool* const filter,
			const bool keep);

};

} /* namespace siena */

#endif /* GENERALTIEITERATOR_H_ */
