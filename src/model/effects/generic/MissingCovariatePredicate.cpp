/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MissingCovariatePredicate.cpp
 *
 * Description: This file contains the implementation of the class
 * MissingCovariatePredicate.
 *****************************************************************************/

#include "MissingCovariatePredicate.h"

using namespace std;

namespace siena
{

/**
 * Constructs a new predicate.
 */
MissingCovariatePredicate::MissingCovariatePredicate(string covariateName) :
	CovariatePredicate(covariateName)
{
}


/**
 * Returns if this predicate holds for the given alter. It is assumed
 * that the predicate has been initialized before and pre-processed with
 * respect to a certain ego.
 */
bool MissingCovariatePredicate::value(int alter)
{
	return this->missing(this->ego()) || this->missing(alter);
}

}
