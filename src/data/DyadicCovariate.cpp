/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: DyadicCovariate.cpp
 * 
 * Description: This file contains the implementation of the
 * DyadicCovariate class.
 *****************************************************************************/

#include "DyadicCovariate.h"

namespace siena
{

/**
 * Constructs a dyadic covariate between a pair of actor sets.
 * @param[in] name the name of the covariate
 * @param[in] pFirstActorSet one of the involved actor sets
 * @param[in] pSecondActorSet the other of the involved actor sets
 */
DyadicCovariate::DyadicCovariate(std::string name,
	const ActorSet * pFirstActorSet,
	const ActorSet * pSecondActorSet) : NamedObject(name)
{
	this->lpFirstActorSet = pFirstActorSet;
	this->lpSecondActorSet = pSecondActorSet;
	this->lmean = 0;
}


/**
 * Returns the first of the actor sets underlying this dyadic covariate.
 */
const ActorSet * DyadicCovariate::pFirstActorSet() const
{
	return this->lpFirstActorSet;
}


/**
 * Returns the second of the actor sets underlying this dyadic covariate.
 */
const ActorSet * DyadicCovariate::pSecondActorSet() const
{
	return this->lpSecondActorSet;
}


/**
 * Stores the average value of this covariate.
 */
void DyadicCovariate::mean(double value)
{
	this->lmean = value;
}

}
