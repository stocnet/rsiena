/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NamedObject.h
 *
 * Description: This file contains the definition of the
 * NamedObject class.
 *****************************************************************************/

#ifndef NAMEDOBJECT_H_
#define NAMEDOBJECT_H_

#include <string>

namespace siena
{

/**
 * This class can be used as a base class for objects having names.
 */
class NamedObject
{
public:
	NamedObject(std::string name);

	std::string name() const;

private:
	// The name of this object
	std::string lname {};
};

}

#endif /*NAMEDOBJECT_H_*/
