/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ChangingDyadicCovariate.h
 *
 * Description: This file contains the definition of the
 * ChangingDyadicCovariate class.
 *****************************************************************************/

#ifndef CHANGINGDYADICCOVARIATE_H_
#define CHANGINGDYADICCOVARIATE_H_

#include <map>
#include <set>
#include "DyadicCovariate.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class DyadicCovariateValueIterator;


// ----------------------------------------------------------------------------
// Section: ChangingDyadicCovariate definition
// ----------------------------------------------------------------------------

/**
 * This class defines a dyadic covariate that changes over time.
 */
class ChangingDyadicCovariate : public DyadicCovariate
{
public:
	ChangingDyadicCovariate(std::string name,
		const ActorSet * pFirstActorSet,
		const ActorSet * pSecondActorSet,
		int observationCount);
	virtual ~ChangingDyadicCovariate();

	double value(int i, int j, int observation) const;
	void value(int i, int j, int observation, double value);
	bool missing(int i, int j, int observation) const;
	void missing(int i, int j, int observation, bool flag);
	DyadicCovariateValueIterator rowValues(int i, int observation,
		bool excludeMissings) const;
	DyadicCovariateValueIterator columnValues(int j, int observation,
		bool excludeMissings) const;
	const std::map <int, double> & rRowValues(int i, int observation) const;

private:
	// A row based representation of non-zero values of the covariate.
	// The value at observation k for a pair (i,j) is stored in
	// lpRowValues[k][i][j].

	std::map<int, double> ** lpRowValues;

	// A column based representation of non-zero values of the covariate.
	// The value at observation k for a pair (i,j) is stored in
	// lpColumnValues[k][j][i].

	std::map<int, double> ** lpColumnValues;

	// A row based representation of missing values. Actor j belongs to
	// lpRowMissings[k][i] if and only if the covariate value for the pair
	// (i,j) is missing at observation k.

	std::set<int> ** lpRowMissings;

	// A column based representation of missing values. Actor i belongs to
	// lpColumnMissings[k][j] if and only if the covariate value for the pair
	// (i,j) is missing at observation k.

	std::set<int> ** lpColumnMissings;

	// The number of observations
	int lobservationCount;

	// empty set
	std::set <int> * lpEmptySet;
};

}

#endif /*CHANGINGDYADICCOVARIATE_H_*/
