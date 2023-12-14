/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DifferentCovariateInStarFunction.h
 *
 * Description: This file contains the definition of the
 * DifferentCovariateInStarFunction class.
 *****************************************************************************/

#ifndef DIFFERENTCOVARIATEINSTARFUNCTION_H_
#define DIFFERENTCOVARIATEINSTARFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"

namespace siena
{

class DifferentCovariateInStarFunction: public CovariateNetworkAlterFunction
{
public:
	DifferentCovariateInStarFunction(std::string networkName,
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

#endif /* DIFFERENTCOVARIATEINSTARFUNCTION_H_ */
