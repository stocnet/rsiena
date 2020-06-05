/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Utils.cpp
 *
 * Description: Defines some utilities for general use.
 *****************************************************************************/

#include "Utils.h"
#include <Rmath.h>
#include <sstream>

namespace siena
{

/**
 * Returns a string representing the given integer.
 */
std::string toString(int i)
{
	std::stringstream out;
	out << i;
	return out.str();
}


/**
 * An identity function simply returning its argument.
 */
double identity(int x)
{
	return x;
}


/**
 * Returns the inverse of (<i>x</i> + 1).
 */
double invertor(int x)
{
	return 1.0 / (x + 1);
}

/**
 * Returns the logarithm of (<i>x</i> + 1).
 */
double logarithmer(int x)
{
	return log((double) (x + 1));
}


/**
 * Creates an exception signaling about the use of an invalid iterator.
 */
InvalidIteratorException::InvalidIteratorException() :
	std::logic_error("An attempt to use an invalid iterator")
{
}

}
