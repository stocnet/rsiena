/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ExogenousEvent.h
 * 
 * Description: This file contains the definition of the
 * ExogenousEvent class.
 *****************************************************************************/

#ifndef EXOGENOUSEVENT_H_
#define EXOGENOUSEVENT_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enumerations
// ----------------------------------------------------------------------------

/**
 * This enumeration defines possible types of composition change.
 */
enum EventType { JOINING, LEAVING };


// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ActorSet;


// ----------------------------------------------------------------------------
// Section: ExogenousEvent class
// ----------------------------------------------------------------------------

/**
 * This class encapsulates the data of a single exogenous event of composition
 * change.
 */
class ExogenousEvent
{
public:
	ExogenousEvent(const ActorSet * pActorSet,
		int actor,
		double time,
		EventType type);
	
	const ActorSet * pActorSet() const;
	int actor() const;
	double time() const;
	EventType type() const;
	
private:
	// The actor set under consideration
	const ActorSet * lpActorSet;
	
	// The actor joining or leaving
	int lactor {};
	
	// The time of composition change (in [0,1])
	double ltime {};
	
	// The type of the change
	EventType ltype {};
};

}

#endif /*EXOGENOUSEVENT_H_*/
