/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantDyadicCovariate.cpp
 *
 * Description: This file contains the implementation of the
 * ConstantDyadicCovariate class.
 *****************************************************************************/

#include "ConstantDyadicCovariate.h"
#include "data/ActorSet.h"
#include "data/DyadicCovariateValueIterator.h"

using namespace std;

namespace siena
{

/**
 * Creates a constant dyadic covariate between a pair of actor sets.
 * @param[in] name the name of the covariate
 * @param[in] pFirstActorSet one of the involved actor sets
 * @param[in] pSecondActorSet the other of the involved actor sets
 */
ConstantDyadicCovariate::ConstantDyadicCovariate(std::string name,
	const ActorSet * pFirstActorSet,
	const ActorSet * pSecondActorSet) :
		DyadicCovariate(name, pFirstActorSet, pSecondActorSet)
{
	this->lpRowValues = new map<int, double>[pFirstActorSet->n()];
	this->lpRowMissings = new set<int>[pFirstActorSet->n()];
	this->lpColumnValues = new map<int, double>[pSecondActorSet->n()];
	this->lpColumnMissings = new set<int>[pSecondActorSet->n()];
}


/**
 * Deallocates this covariate object.
 */
ConstantDyadicCovariate::~ConstantDyadicCovariate()
{
	delete[] this->lpRowValues;
	delete[] this->lpRowMissings;
	delete[] this->lpColumnValues;
	delete[] this->lpColumnMissings;
	this->lpRowValues = 0;
	this->lpRowMissings = 0;
	this->lpColumnValues = 0;
	this->lpColumnMissings = 0;
}


/**
 * Stores the value for the given pair of actors.
 * @param[in] i the first actor of the pair
 * @param[in] j the second actor of the pair
 * @param[in] value the value to be stored
 */
void ConstantDyadicCovariate::value(int i, int j, double value)
{
	if (value)
	{
		this->lpRowValues[i][j] = value;
		this->lpColumnValues[j][i] = value;
	}
	else
	{
		this->lpRowValues[i].erase(j);
		this->lpColumnValues[j].erase(i);
	}
}


/**
 * Returns the value of the covariate for the given pair of actors.
 */
double ConstantDyadicCovariate::value(int i, int j) const
{
	map<int, double>::const_iterator iter = this->lpRowValues[i].find(j);
	double value = 0;

	if (iter != this->lpRowValues[i].end())
	{
		value = iter->second;
	}

	return value;
}


/**
 * Stores if the value for the given pair of actors is missing.
 * @param[in] i the first actor of the pair
 * @param[in] j the second actor of the pair
 * @param[in] flag indicates if the value is missing
 */
void ConstantDyadicCovariate::missing(int i, int j, bool flag)
{
	if (flag)
	{
		this->lpRowMissings[i].insert(j);
		this->lpColumnMissings[j].insert(i);
	}
	else
	{
		this->lpRowMissings[i].erase(j);
		this->lpColumnMissings[j].erase(i);
	}
}


/**
 * Returns if the value for the given pair of actors is missing.
 */
bool ConstantDyadicCovariate::missing(int i, int j) const
{
	return this->lpRowMissings[i].find(j) != this->lpRowMissings[i].end();
}


/**
 * Returns an iterator over non-zero non-missing values of the given row.
 */
DyadicCovariateValueIterator ConstantDyadicCovariate::rowValues(int i)
	const
{
	return DyadicCovariateValueIterator(this->lpRowValues[i],
		this->lpRowMissings[i]);
}


/**
 * Returns an iterator over non-zero non-missing values of the given column.
 */
DyadicCovariateValueIterator ConstantDyadicCovariate::columnValues(int j)
	const
{
	return DyadicCovariateValueIterator(this->lpColumnValues[j],
		this->lpColumnMissings[j]);
}

/**
 * Returns a map containing over non-zero non-missing values of the given row.
 */
const map <int, double> & ConstantDyadicCovariate::rRowValues(int i)
	const
{
	return this->lpRowValues[i];
}
}
