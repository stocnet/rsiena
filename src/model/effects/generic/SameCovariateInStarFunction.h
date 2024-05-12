/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateInStarFunction.h
 *
 * Description: This file contains the definition of the
 * SameCovariateInStarFunction class.
 *****************************************************************************/

#ifndef SAMECOVARIATEINSTARFUNCTION_H_
#define SAMECOVARIATEINSTARFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"

namespace siena
{

class SameCovariateInStarFunction: public CovariateNetworkAlterFunction
{
public:
	SameCovariateInStarFunction(std::string networkName,
		std::string covariateName, bool excludeMissing);
		
	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
};

}

#endif /* SAMECOVARIATEINSTARFUNCTION_H_ */
