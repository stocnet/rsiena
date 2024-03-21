/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ChangingCovariate.cpp
 * 
 * Description: This file contains the implementation of the
 * ChangingCovariate class.
 *****************************************************************************/

#include <limits>
#include "ChangingCovariate.h"
#include "data/ActorSet.h"

using namespace std;

namespace siena
{

/**
 * Creates a new changing covariate for the given set of actors and the
 * given number of observations. Initially all values are set to 0.
 * @param[in] name the name of the covariate
 */
ChangingCovariate::ChangingCovariate(std::string name,
	const ActorSet * pActorSet,
	int observationCount) : Covariate(name, pActorSet)
{
	this->lvalues = new double * [pActorSet->n()];
	this->lmissing = new bool * [pActorSet->n()];
	
	for (int i = 0; i < pActorSet->n(); i++)
	{
		this->lvalues[i] = new double[observationCount];
		this->lmissing[i] = new bool[observationCount];
		
		for (int j = 0; j < observationCount; j++)
		{
			this->lvalues[i][j] = 0;
			this->lmissing[i][j] = false;
		}
	}
	this->lmin = numeric_limits<double>::max();
	this->lmax = numeric_limits<double>::min();
}


/**
 * Deallocates this changing covariate object.
 */
ChangingCovariate::~ChangingCovariate()
{
	for (int i = 0; i < this->pActorSet()->n(); i++)
	{
		delete[] this->lvalues[i];
		delete[] this->lmissing[i];
	}
	
	delete[] this->lvalues;
	delete[] this->lmissing;
	this->lvalues = 0;
	this->lmissing = 0;
}


/**
 * Returns the value of this covariate for actor <i>i</i> at the given
 * observation.
 */
double ChangingCovariate::value(int i, int observation) const
{
	return this->lvalues[i][observation];
}


/**
 * Stores the value of this covariate for actor <i>i</i> at the given
 * observations.
 */
void ChangingCovariate::value(int i, int observation, double value)
{
	this->lvalues[i][observation] = value;
	this->lmin = std::min(this->lmin, value);
	this->lmax = std::max(this->lmax, value);
}


/**
 * Returns if the value of the covariate is missing for the given
 * actor at the specified observation.
 */
bool ChangingCovariate::missing(int i, int observation) const
{
	return this->lmissing[i][observation];
}


/**
 * Stores if the value of the covariate is missing for the given
 * actor at the specified observation.
 */
void ChangingCovariate::missing(int i,
	int observation,
	bool missing)
{
	this->lmissing[i][observation] = missing;
}

/**
 * Returns the smallest observed value.
 */
double ChangingCovariate::min() const
{
	return this->lmin;
}

/**
 * Returns the largest observed value.
 */
double ChangingCovariate::max() const
{
	return this->lmax;
}
}
