/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ChangingDyadicCovariate.cpp
 *
 * Description: This file contains the implementation of the
 * ChangingDyadicCovariate class.
 *****************************************************************************/

#include <R_ext/Print.h>
#include "ChangingDyadicCovariate.h"
#include "data/ActorSet.h"
#include "data/DyadicCovariateValueIterator.h"

using namespace std;

namespace siena
{

/**
 * Creates a changing dyadic covariate between a pair of actor sets.
 * @param[in] name the name of the covariate
 * @param[in] pFirstActorSet one of the involved actor sets
 * @param[in] pSecondActorSet the other of the involved actor sets
 * @param[in] observationCount the number of observations of the covariate
 */
ChangingDyadicCovariate::ChangingDyadicCovariate(string name,
	const ActorSet * pFirstActorSet,
	const ActorSet * pSecondActorSet,
	int observationCount) :
		DyadicCovariate(name, pFirstActorSet, pSecondActorSet)
{
	this->lobservationCount = observationCount;
	this->lpRowValues = new map<int, double> * [observationCount];
	this->lpColumnValues = new map<int, double> * [observationCount];
	this->lpRowMissings = new set<int> * [observationCount];
	this->lpColumnMissings = new set<int> * [observationCount];

	for (int k = 0; k < observationCount; k++)
	{
		this->lpRowValues[k] = new map<int, double>[pFirstActorSet->n()];
		this->lpColumnValues[k] = new map<int, double>[pSecondActorSet->n()];
		this->lpRowMissings[k] = new set<int>[pFirstActorSet->n()];
		this->lpColumnMissings[k] = new set<int>[pSecondActorSet->n()];
	}
	this->lpEmptySet = new set<int>;
}


/**
 * Deallocates this covariate object.
 */
ChangingDyadicCovariate::~ChangingDyadicCovariate()
{
	for (int k = 0; k < this->lobservationCount; k++)
	{
		delete[] this->lpRowValues[k];
		delete[] this->lpColumnValues[k];
		delete[] this->lpRowMissings[k];
		delete[] this->lpColumnMissings[k];
	}

	delete[] this->lpRowValues;
	delete[] this->lpColumnValues;
	delete[] this->lpRowMissings;
	delete[] this->lpColumnMissings;
	this->lpRowValues = 0;
	this->lpColumnValues = 0;
	this->lpRowMissings = 0;
	this->lpColumnMissings = 0;
	delete this->lpEmptySet;
	this->lpEmptySet = 0;
}


/**
 * Stores the value for the given pair of actors at the given observation.
 * @param[in] i the first actor of the pair
 * @param[in] j the second actor of the pair
 * @param[in] observation the number of the observation
 * @param[in] value the value to be stored
 */
void ChangingDyadicCovariate::value(int i,
	int j,
	int observation,
	double value)
{
	if (value)
	{
		this->lpRowValues[observation][i][j] = value;
		this->lpColumnValues[observation][j][i] = value;
	}
	else
	{
		this->lpRowValues[observation][i].erase(j);
		this->lpColumnValues[observation][j].erase(i);
	}
}


/**
 * Returns the value of the covariate for the given pair of actors at the
 * given observation.
 */
double ChangingDyadicCovariate::value(int i, int j, int observation) const
{
	map<int, double>::const_iterator iter =
		this->lpRowValues[observation][i].find(j);
	double value = 0;

	if (iter != this->lpRowValues[observation][i].end())
	{
		value = iter->second;
	}

	return value;
}


/**
 * Stores if the value for the given pair of actors is missing
 * at the given observation.
 * @param[in] i the first actor of the pair
 * @param[in] j the second actor of the pair
 * @param[in] observation the number of the observation
 * @param[in] flag indicates if the value is missing
 */
void ChangingDyadicCovariate::missing(int i,
	int j,
	int observation,
	bool flag)
{
	if (flag)
	{
		this->lpRowMissings[observation][i].insert(j);
		this->lpColumnMissings[observation][j].insert(i);
	}
	else
	{
		this->lpRowMissings[observation][i].erase(j);
		this->lpColumnMissings[observation][j].erase(i);
	}
}


/**
 * Returns if the value is missing for the given pair of actors at the
 * given observation.
 */
bool ChangingDyadicCovariate::missing(int i, int j, int observation) const
{
	return this->lpRowMissings[observation][i].find(j) !=
		this->lpRowMissings[observation][i].end();
}


/**
 * Returns an iterator over non-zero values of the given row
 * at the given observation.
 * @param[in] excludeMissings indicates if missing values should be
 * excluded from the iteration
 */
DyadicCovariateValueIterator ChangingDyadicCovariate::rowValues(int i,
	int observation,
	bool excludeMissings) const
{
	set<int> * excludedActorSet = this->lpEmptySet;

    if (excludeMissings)
    {
        excludedActorSet = &(this->lpRowMissings[observation][i]);
    }

    return DyadicCovariateValueIterator(this->lpRowValues[observation][i],
        *excludedActorSet);
}


/**
 * Returns an iterator over non-zero values of the given column
 * at the given observation.
 * @param[in] excludeMissings indicates if missing values should be
 * excluded from the iteration
 */
DyadicCovariateValueIterator ChangingDyadicCovariate::columnValues(int j,
	int observation,
	bool excludeMissings) const
{
	set<int> * excludedActorSet = this->lpEmptySet;
//	Rprintf("%d exclude\n", excludeMissings);

    if (excludeMissings)
    {
		//	Rprintf("exclude\n");
        excludedActorSet = &(this->lpColumnMissings[observation][j]);
    }

    return DyadicCovariateValueIterator(this->lpColumnValues[observation][j],
        *excludedActorSet);
}
/**
 * Returns a map containing over non-zero non-missing values of the given row
 * for the given observation
 */
const map <int, double> & ChangingDyadicCovariate::rRowValues(int i,
	int observation) const
{
	return this->lpRowValues[observation][i];
}

}
