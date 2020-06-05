/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: NamedObject.cpp
 * 
 * Description: This file contains the implementation of the
 * NamedObject class.
 *****************************************************************************/

#include "NamedObject.h"

namespace siena
{

/**
 * Constructs an object with the given name.
 */
NamedObject::NamedObject(std::string name)
{
	this->lname = name;
}


/**
 * Returns the name of this object.
 */
std::string NamedObject::name() const
{
	return this->lname;
}

}
