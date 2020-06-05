/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateValueIterator.h
 *
 * Description: This file contains the definition of the
 * DyadicCovariateValueIterator class.
 *****************************************************************************/

#ifndef DYADICCOVARIATEVALUEITERATOR_H_
#define DYADICCOVARIATEVALUEITERATOR_H_

#include <map>
#include <set>

namespace siena
{

/**
 * Defines an iterator over non-zero values w_{ij} of a dyadic covariate
 * for a specific actor i.
 */
class DyadicCovariateValueIterator
{
	// These classes needs access to the private constructor.

	friend class ConstantDyadicCovariate;
	friend class ChangingDyadicCovariate;

public:
	int actor() const;
	double value() const;
	bool valid() const;
	void next();

private:
	DyadicCovariateValueIterator(std::map<int, double> & rValues,
		std::set<int> & rMissings);
	void skipMissings();

	// Points to the current element in the map of values
	std::map<int, double>::const_iterator lcurrent;

	// Points to the end of the map of values
	std::map<int, double>::const_iterator lend;

	// Points to the current element in the set of missing values
	std::set<int>::const_iterator lmissingCurrent;

	// Points to the end of the set of missing values
	std::set<int>::const_iterator lmissingEnd;
};

}

#endif /*DYADICCOVARIATEVALUEITERATOR_H_*/
