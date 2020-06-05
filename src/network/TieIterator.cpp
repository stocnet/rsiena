/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: TieIterator.cpp
 * 
 * Description: This file contains the implementation of the
 * TieIterator class.
 *****************************************************************************/

#include "TieIterator.h"
#include "Network.h"
#include "../utils/Utils.h"


namespace siena
{

/**
 * Creates an iterator over all ties of the given network.
 */
TieIterator::TieIterator(const Network * pNetwork)
{
	this->lpNetwork = pNetwork;
	
	if (this->lpNetwork->n() == 0)
	{
		// No actors means no ties. There is nothing to iterate over.
		this->lvalid = false;
	}
	else
	{
		this->lvalid = true;
		this->lcurrentActor = 0;
		this->literator = this->lpNetwork->outTies(0);
		
		// Find the first actor with some outgoing ties. If there is no
		// such an actor, the iterator will be invalidated.
		
		this->skipInvalidIterators();
	}
}

//
// If the outgoing tie iterator of the current actor has become invalid,
// the current actor is incremented until we reach an actor with some
// outgoing ties or run out of actors, in which case this iterator becomes
// invalid.
//
void TieIterator::skipInvalidIterators()
{
	while (this->lvalid && !this->literator.valid())
	{
		this->lcurrentActor++;
		
		if (this->lcurrentActor == this->lpNetwork->n())
		{
			// No more actors. We have run out of ties.
			this->lvalid = false;
		}
		else
		{
			this->literator = this->lpNetwork->outTies(this->lcurrentActor);
		}
	}
}


/**
 * Returns the sender of the current tie.
 */
int TieIterator::ego() const
{
	this->checkValidity();
	return this->lcurrentActor;
}


/**
 * Returns the receiver of the current tie.
 */
int TieIterator::alter() const
{
	this->checkValidity();
	return this->literator.actor();
}


/**
 * Returns the value of the current tie.
 */
int TieIterator::value() const
{
	this->checkValidity();
	return this->literator.value();
}


/**
 * Indicates if this iterator still points to a valid tie.
 */
bool TieIterator::valid() const
{
	return this->lvalid;
}


/**
 * Advances the iterator to the next tie.
 */
void TieIterator::next()
{
	this->checkValidity();
	this->literator.next();
	this->skipInvalidIterators();
}


/**
 * Throws an InvalidIteratorException if this iterator is invalid.
 */
void TieIterator::checkValidity() const
{
	if (!this->valid())
	{
		throw InvalidIteratorException();
	}
}

}
