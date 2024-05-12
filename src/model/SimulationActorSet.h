/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimulationActorSet.h
 *
 * Description: This file contains the definition of the
 * SimulationActorSet class.
 *****************************************************************************/

#ifndef SIMULATIONACTORSET_H_
#define SIMULATIONACTORSET_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ActorSet;


// ----------------------------------------------------------------------------
// Section: Class SimulationActorSet
// ----------------------------------------------------------------------------

/**
 * A wrapper of actor sets for simulation purposes.
 */
class SimulationActorSet
{
public:
	SimulationActorSet(const ActorSet * pActorSet);
	virtual ~SimulationActorSet();

	inline const ActorSet * pOriginalActorSet() const;
	inline int n() const;
	inline bool active(int i) const;
	inline int activeActorCount() const;

	void active(int i, bool flag);

private:
	// The real actor set represented by this simulation actor set
	const ActorSet * lpActorSet;

	// The number of actors in this set
	int ln {};

	// An indicator per actor if the actor is currently active
	bool * lactive {};

	// The number of currently active actors
	int lactiveActorCount {};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the original actor set this object represents.
 */
const ActorSet * SimulationActorSet::pOriginalActorSet() const
{
	return this->lpActorSet;
}


/**
 * Returns the number of actors in this actor set.
 */
int SimulationActorSet::n() const
{
	return this->ln;
}


/**
 * Returns if the given actor is currently active.
 */
bool SimulationActorSet::active(int i) const
{
	return this->lactive[i];
}


/**
 * Returns the number of currently active actors.
 */
int SimulationActorSet::activeActorCount() const
{
	return this->lactiveActorCount;
}

}

#endif /* SIMULATIONACTORSET_H_ */
