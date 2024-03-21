/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleEqualCovariateFunction.h
 *
 * Description: This file contains the definition of the
 * DoubleEqualCovariateFunction class.
 *****************************************************************************/

#ifndef DOUBLEEQUALCOVARIATEFUNCTION_H_
#define DOUBLEEQUALCOVARIATEFUNCTION_H_

#include "DoubleCovariateFunction.h"

namespace siena
{

/**
 * Defines a predicate that holds if the ego and the alter have the same
 * covariate values.
 */
class DoubleEqualCovariateFunction: public DoubleCovariateFunction
{
public:
	DoubleEqualCovariateFunction(std::string covariateName1,
						std::string covariateName2, bool excludeMissing);
						
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
};

}

#endif /* DOUBLEEQUALCOVARIATEFUNCTION_H_ */
