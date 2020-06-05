/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ActorSet.cpp
 *
 * Description: This file contains the implementation of the
 * ActorSet class.
 *****************************************************************************/

#include "ActorSet.h"

namespace siena
{

/**
 * Creates a new set of actors of the given size.
 */
ActorSet::ActorSet(std::string name, int n) : NamedObject(name)
{
	this->ln = n;
}

/**
 * Deallocates this set of actors.
 */
ActorSet::~ActorSet()
{
	//Rprintf("delete actor Set\n");
}
}
