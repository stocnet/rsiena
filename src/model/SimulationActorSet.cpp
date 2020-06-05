/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimulationActorSet.cpp
 *
 * Description: This file contains the implementation of the class
 * SimulationActorSet.
 *****************************************************************************/

#include "SimulationActorSet.h"
#include "data/ActorSet.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pActorSet the real actor set that this object is to represent
 */
SimulationActorSet::SimulationActorSet(const ActorSet * pActorSet)
{
	this->lpActorSet = pActorSet;
	this->ln = pActorSet->n();

	this->lactive = new bool[this->ln];

	for (int i = 0; i < this->ln; i++)
	{
		this->lactive[i] = false;
	}

	this->lactiveActorCount = 0;
}


/**
 * Destructor.
 */
SimulationActorSet::~SimulationActorSet()
{
	delete[] this->lactive;
	this->lactive = 0;
}


/**
 * Sets the activity flag of the given actor.
 */
void SimulationActorSet::active(int i, bool flag)
{
	if (this->lactive[i] != flag)
	{
		this->lactive[i] = flag;

		if (flag)
		{
			// The actor is becoming active
			this->lactiveActorCount++;
		}
		else
		{
			// The actor is going inactive
			this->lactiveActorCount--;
		}
	}
}

}
