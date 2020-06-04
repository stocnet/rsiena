/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EqualCovariatePredicate.h
 *
 * Description: This file contains the definition of the
 * EqualCovariatePredicate class.
 *****************************************************************************/

#ifndef EQUALCOVARIATEPREDICATE_H_
#define EQUALCOVARIATEPREDICATE_H_

#include "CovariatePredicate.h"

namespace siena
{

/**
 * Defines a predicate that holds if the ego and the alter have the same
 * covariate values.
 */
class EqualCovariatePredicate: public CovariatePredicate
{
public:
	EqualCovariatePredicate(std::string covariateName);

	virtual bool value(int alter);
};

}

#endif /* EQUALCOVARIATEPREDICATE_H_ */
