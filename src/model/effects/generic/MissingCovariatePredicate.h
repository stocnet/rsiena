/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MissingCovariatePredicate.h
 *
 * Description: This file contains the definition of the
 * MissingCovariatePredicate class.
 *****************************************************************************/

#ifndef MISSINGCOVARIATEPREDICATE_H_
#define MISSINGCOVARIATEPREDICATE_H_

#include "CovariatePredicate.h"

namespace siena
{

/**
 * Defines a predicate that holds if the covariate value is missing either for
 * the ego or the alter.
 */
class MissingCovariatePredicate: public CovariatePredicate
{
public:
	MissingCovariatePredicate(std::string covariateName);

	virtual bool value(int alter);
};

}

#endif /* MISSINGCOVARIATEPREDICATE_H_ */
