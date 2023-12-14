/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DifferentCovariateOutStarFunction.h
 *
 * Description: This file contains the definition of the
 * DifferentCovariateOutStarFunction class.
 *****************************************************************************/

#ifndef DIFFERENTCOVARIATEOUTSTARFUNCTION_H_
#define DIFFERENTCOVARIATEOUTSTARFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"

namespace siena
{

class DifferentCovariateOutStarFunction: public CovariateNetworkAlterFunction
{
public:
	DifferentCovariateOutStarFunction(std::string networkName,
		std::string covariateName, bool excludeMissing, bool bothDifferent);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

		virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
	bool lnotBothDifferent {};
};

}

#endif /* DIFFERENTCOVARIATEOUTSTARFUNCTION_H_ */
