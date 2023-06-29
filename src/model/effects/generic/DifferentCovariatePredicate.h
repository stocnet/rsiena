/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DifferentCovariatePredicate.h
 *
 * Description: This file contains the definition of the
 * DifferentCovariatePredicate class.
 *****************************************************************************/

#ifndef DIFFERENTCOVARIATEPREDICATE_H_
#define DIFFERENTCOVARIATEPREDICATE_H_

#include "CovariatePredicate.h"

namespace siena
{

/**
 * Defines a predicate that holds if the ego and the alter have the same
 * covariate values.
 */
class DifferentCovariatePredicate: public CovariatePredicate
{
public:
	DifferentCovariatePredicate(std::string covariateName);

	virtual bool value(int alter);
};

}

#endif /* DIFFERENTCOVARIATEPREDICATE_H_ */
