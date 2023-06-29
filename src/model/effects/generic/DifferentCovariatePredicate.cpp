/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DifferentCovariatePredicate.cpp
 *
 * Description: This file contains the implementation of the class
 * DifferentCovariatePredicate.
 *****************************************************************************/

#include <cmath>
#include "DifferentCovariatePredicate.h"
#include "utils/Utils.h"

using namespace std;

namespace siena
{

/**
 * Constructs a new predicate.
 */
DifferentCovariatePredicate::DifferentCovariatePredicate(string covariateName) :
	CovariatePredicate(covariateName)
{
}

/**
 * Returns if this predicate holds for the given alter. It is assumed
 * that the predicate has been initialized before and pre-processed with
 * respect to a certain ego.
 */
bool DifferentCovariatePredicate::value(int alter) 
{
	return fabs(this->covariateValue(this->ego()) -
		this->covariateValue(alter)) >= EPSILON;
}

}
