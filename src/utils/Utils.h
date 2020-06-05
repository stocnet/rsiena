/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Utils.h
 *
 * Description: Defines some utilities for general use.
 *****************************************************************************/

#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <stdexcept>
#include <list>
#include <vector>
#include <map>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Macros
// ----------------------------------------------------------------------------

#define EPSILON 1e-6

// ----------------------------------------------------------------------------
// Section: Methods
// ----------------------------------------------------------------------------

std::string toString(int i);
double identity(int x);
double invertor(int x);
double logarithmer(int x);

/**
 * This method tests if the given element belongs to the given container.
 */
template<class C, class E>
bool contains(const C & rContainer, const E element)
{
	return find(rContainer.begin(), rContainer.end(), element) !=
		rContainer.end();
}


/**
 * Deletes each element of the given list and empties the list itself.
 */
template<class T>
void deallocateList(std::list<T *> & rList)
{
	while (!rList.empty())
	{
		T * pItem = rList.front();
		rList.pop_front();
		delete pItem;
	}

	rList.clear();
}


/**
 * Deletes each element of the given vector and empties the vector itself.
 */
template<class T>
void deallocateVector(std::vector<T *> & rVector)
{
	for (unsigned i = 0; i < rVector.size(); i++)
	{
		delete rVector[i];
	}

	rVector.clear();
}


/**
 * Removes all elements from the given map. The keys and/or values are
 * deleted on request.
 */
template<class K, class V>
void clearMap(std::map<K *, V *> & rMap, bool deleteKeys, bool deleteValues)
{
	while (!rMap.empty())
	{
		K * pKey = rMap.begin()->first;
		V * pValue = rMap.begin()->second;

		rMap.erase(rMap.begin());

		if (deleteKeys)
		{
			delete pKey;
		}

		if (deleteValues)
		{
			delete pValue;
		}
	}
}


// ----------------------------------------------------------------------------
// Section: Exception classes
// ----------------------------------------------------------------------------

/**
 * This class defines an exception signaling about the use of an invalid
 * iterator.
 */
class InvalidIteratorException : public std::logic_error
{
public:
	InvalidIteratorException();
};

}

#endif /*UTILS_H_*/
