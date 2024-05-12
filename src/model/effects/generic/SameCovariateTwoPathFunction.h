/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateTwoPathFunction.h
 *
 * Description: This file contains the definition of the
 * SameCovariateTwoPathFunction class.
 *****************************************************************************/

#ifndef SAMECOVARIATETWOPATHFUNCTION_H_
#define SAMECOVARIATETWOPATHFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"

namespace siena
{

class SameCovariateTwoPathFunction: public CovariateNetworkAlterFunction
{
public:
	SameCovariateTwoPathFunction(std::string networkName,
		std::string covariateName, bool same, bool excludeMissing);
		
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lsame {};
	bool lexcludeMissing {};
};

}

#endif /* SAMECOVARIATETWOPATHFUNCTION_H_ */
