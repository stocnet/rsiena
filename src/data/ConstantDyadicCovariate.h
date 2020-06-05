/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantDyadicCovariate.h
 *
 * Description: This file contains the definition of the
 * ConstantDyadicCovariate class.
 *****************************************************************************/

#ifndef CONSTANTDYADICCOVARIATE_H_
#define CONSTANTDYADICCOVARIATE_H_

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
// Section: ConstantDyadicCovariate class
// ----------------------------------------------------------------------------

/**
 * This class defines a constant dyadic covariate.
 */
class ConstantDyadicCovariate : public DyadicCovariate
{
public:
	ConstantDyadicCovariate(std::string name,
		const ActorSet * pFirstActorSet,
		const ActorSet * pSecondActorSet);
	virtual ~ConstantDyadicCovariate();

	double value(int i, int j) const;
	void value(int i, int j, double value);
	bool missing(int i, int j) const;
	void missing(int i, int j, bool flag);
	DyadicCovariateValueIterator rowValues(int i) const;
	DyadicCovariateValueIterator columnValues(int j) const;
	const std::map <int, double> & rRowValues(int i) const;

private:
	// A row based representation of non-zero values of the covariate.
	// The value for a pair (i,j) is stored in lpRowValues[i][j].

	std::map<int, double> * lpRowValues;

	// A column based representation of non-zero values of the covariate.
	// The value for a pair (i,j) is stored in lpColumnValues[j][i].

	std::map<int, double> * lpColumnValues;

	// A row based representation of missing values. Actor j belongs to
	// lpRowMissings[i] if and only if the covariate value for the pair
	// (i,j) is missing.

	std::set<int> * lpRowMissings;

	// A column based representation of missing values. Actor i belongs to
	// lpColumnMissings[j] if and only if the covariate value for the pair
	// (i,j) is missing.

	std::set<int> * lpColumnMissings;
};

}

#endif /*CONSTANTDYADICCOVARIATE_H_*/
