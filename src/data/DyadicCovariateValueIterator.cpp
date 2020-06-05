/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateValueIterator.cpp
 *
 * Description: This file contains the implementation of the class
 * DyadicCovariateValueIterator.
 *****************************************************************************/

#include "DyadicCovariateValueIterator.h"
#include "utils/Utils.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] rValues the map storing the non-zero values w_{ij} for each j
 * @param[in] rMissings a set of actors for which the values w_{ij} are missing
 */
DyadicCovariateValueIterator::DyadicCovariateValueIterator(
	map<int, double> & rValues,
	set<int> & rMissings)
{
	this->lcurrent = rValues.begin();
	this->lend = rValues.end();
	this->lmissingCurrent = rMissings.begin();
	this->lmissingEnd = rMissings.end();

	this->skipMissings();
}


/**
 * Returns the current actor j with a non-zero non-missing value
 * w_{ij} of the covariate.
 */
int DyadicCovariateValueIterator::actor() const
{
	if (!this->valid())
	{
		throw InvalidIteratorException();
	}

	return this->lcurrent->first;
}


/**
 * Returns the current non-zero non-missing value w_{ij} of the covariate.
 */
double DyadicCovariateValueIterator::value() const
{
	if (!this->valid())
	{
		throw InvalidIteratorException();
	}

	return this->lcurrent->second;
}


/**
 * Indicates if the iterator still points to a valid actor j with a non-zero
 * non-missing covariate value w_{ij}.
 */
bool DyadicCovariateValueIterator::valid() const
{
	return this->lcurrent != this->lend;
}


/**
 * Moves the iterator to the next actor j with a non-zero non-missing covariate
 * value w_{ij}.
 */
void DyadicCovariateValueIterator::next()
{
	this->lcurrent++;
	this->skipMissings();
}


/**
 * Makes sure that lcurrent points to the next non-zero value of the
 * covariate that is not missing.
 */
void DyadicCovariateValueIterator::skipMissings()
{
	while (this->lmissingCurrent != this->lmissingEnd &&
		this->lcurrent != this->lend &&
		(*this->lmissingCurrent) <= this->lcurrent->first)
	{
		if ((*this->lmissingCurrent) == this->lcurrent->first)
		{
			// The current iterator points to a non-zero value,
			// but it is missing, so we skip it.

			this->lcurrent++;
		}

		this->lmissingCurrent++;
	}
}

}
