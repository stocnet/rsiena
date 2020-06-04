/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EqualCovariatePredicate.cpp
 *
 * Description: This file contains the implementation of the class
 * EqualCovariatePredicate.
 *****************************************************************************/

#include <cmath>
#include "EqualCovariatePredicate.h"
#include "utils/Utils.h"

using namespace std;

namespace siena
{

/**
 * Constructs a new predicate.
 */
EqualCovariatePredicate::EqualCovariatePredicate(string covariateName) :
	CovariatePredicate(covariateName)
{
}

/**
 * Returns if this predicate holds for the given alter. It is assumed
 * that the predicate has been initialized before and pre-processed with
 * respect to a certain ego.
 */
bool EqualCovariatePredicate::value(int alter)
{
	return fabs(this->covariateValue(this->ego()) -
		this->covariateValue(alter)) < EPSILON;
}

}
