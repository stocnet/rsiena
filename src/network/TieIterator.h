/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: TieIterator.h
 * 
 * Description: This module defines the class TieIterator for
 * convenient iteration over all ties of a network.
 *****************************************************************************/

#ifndef TIEITERATOR_H_
#define TIEITERATOR_H_

#include "IncidentTieIterator.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;


// ----------------------------------------------------------------------------
// Section: TieIterator class
// ----------------------------------------------------------------------------

/**
 * This class defines an iterator over all ties of a network. The iterator
 * first returns all outgoing ties of actor 0, followed by all outgoing ties of
 * actor 1, and so on.
 */
class TieIterator
{
public:
	TieIterator(const Network * pNetwork);

	int ego() const;
	int alter() const;
	int value() const;
	bool valid() const;
	void next();
	
private:
	void skipInvalidIterators();
	void checkValidity() const;
	
	// The underlying network
	const Network * lpNetwork;
	
	// The actor whose outgoing ties are currently iterated over
	int lcurrentActor;

	// The iterator over outgoing ties of the current actor
	IncidentTieIterator literator;
	
	// Indicates if the iterator still points to an existing tie
	int lvalid {};
};

}

#endif /*TIEITERATOR_H_*/
