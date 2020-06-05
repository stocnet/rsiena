/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ExogenousEvent.cpp
 * 
 * Description: This file contains the implementation of the
 * ExogenousEvent class.
 *****************************************************************************/

#include "ExogenousEvent.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructor area
// ----------------------------------------------------------------------------

/**
 * Creates an exogenous event.
 * @param[in] pActorSet the set of actors affected by this event
 * @param[in] actor the actor joining or leaving
 * @param[in] time the time of composition change
 * @param[in] type the type of the change (joiner of leaver)
 */
ExogenousEvent::ExogenousEvent(const ActorSet * pActorSet,
	int actor,
	double time,
	EventType type)
{
	this->lpActorSet = pActorSet;
	this->lactor = actor;
	this->ltime = time;
	this->ltype = type;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the set of actors affected by this event.
 */
const ActorSet * ExogenousEvent::pActorSet() const
{
	return this->lpActorSet;
}


/**
 * Returns the actor involved in this composition change event.
 */
int ExogenousEvent::actor() const
{
	return this->lactor;
}


/**
 * Returns the time of this composition change. The value is in the range
 * [0,1].
 */
double ExogenousEvent::time() const
{
	return this->ltime;
}


/**
 * Returns the type of this composition change (joiner or leaver).
 */
EventType ExogenousEvent::type() const
{
	return this->ltype;
}

}
