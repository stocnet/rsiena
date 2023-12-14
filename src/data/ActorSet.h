/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ActorSet.h
 *
 * Description: This file contains the definition of the
 * ActorSet class.
 *****************************************************************************/

#ifndef ACTORSET_H_
#define ACTORSET_H_

#include "utils/NamedObject.h"

namespace siena
{

/**
 * This class defines a set of actors of a certain size <i>n</i>. The actors
 * are numbered from 0 to <i>n</i> - 1.
 */
class ActorSet : public NamedObject
{
public:
	ActorSet(std::string name, int n);
	virtual ~ActorSet();

	inline int n() const;

private:
	// The number of actors in this set
	int ln{};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the number of actors in this set.
 */
int ActorSet::n() const
{
	return this->ln;
}

}

#endif /*ACTORSET_H_*/
