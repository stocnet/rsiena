#ifndef UNIONNEIGHBORITERATOR_H_
#define UNIONNEIGHBORITERATOR_H_

#include "IncidentTieIterator.h"
#include "iterators/ITieIterator.h"

namespace siena {

/**
 * This class defines an iterator over actors that are in either of to a pair
 * of incident tie iterators.
 */
class UnionNeighborIterator: ITieIterator {
public:
	UnionNeighborIterator(const IncidentTieIterator& iter1,
			const IncidentTieIterator& iter2);
	UnionNeighborIterator(const UnionNeighborIterator& rhs);
	bool valid() const;
	int actor() const;
	void next();
	UnionNeighborIterator* clone() const;
	// added for GIT version
	virtual void reset();
	int size() const {
	  return -1;
	}
	
private:
	IncidentTieIterator liter1;
	IncidentTieIterator liter2;

	UnionNeighborIterator& operator=(const UnionNeighborIterator&);
};

}

#endif /*UNIONNEIGHBORITERATOR_H_*/
