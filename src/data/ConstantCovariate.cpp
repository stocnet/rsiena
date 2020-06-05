/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ConstantCovariate.cpp
 * 
 * Description: This file contains the implementation of the
 * ConstantCovariate class.
 *****************************************************************************/

#include "ConstantCovariate.h"
#include "data/ActorSet.h"

namespace siena
{

/**
 * Creates a new constant covariate with the given name for the given set
 * of actors.
 */
ConstantCovariate::ConstantCovariate(std::string name,
	const ActorSet * pActorSet) : Covariate(name, pActorSet)
{
	this->lvalues = new double[pActorSet->n()];
	this->lmissing = new bool[pActorSet->n()];
}


/**
 * Deallocates this covariate object.
 */
ConstantCovariate::~ConstantCovariate()
{
	delete[] this->lvalues;
	delete[] this->lmissing;
	this->lvalues = 0;
	this->lmissing = 0;
}


/**
 * Returns the value of this covariate for actor <i>i</i>.
 */
double ConstantCovariate::value(int i) const
{
	return this->lvalues[i];
}


/**
 * Stores the value of this covariate for actor <i>i</i>.
 */
void ConstantCovariate::value(int i, double value)
{
	this->lvalues[i] = value;
}


/**
 * Returns if the value of the covariate is missing for the given actor.
 */
bool ConstantCovariate::missing(int i) const
{
	return this->lmissing[i];
}


/**
 * Stores if the value of the covariate is missing for the given actor.
 */
void ConstantCovariate::missing(int i, bool missing)
{
	this->lmissing[i] = missing;
}

}
