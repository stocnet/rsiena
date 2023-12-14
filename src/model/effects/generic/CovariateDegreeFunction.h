/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDegreeFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDegreeFunction class.
 *****************************************************************************/

#ifndef COVARIATEDEGREEFUNCTION_H_
#define COVARIATEDEGREEFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"

namespace siena
{

class CovariateDegreeFunction: public CovariateNetworkAlterFunction
{
public:
	CovariateDegreeFunction(std::string networkName,
		std::string covariateName, bool excludeMissing,
						bool incoming, bool forEgo, bool sqrtVersion);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

		virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
	bool lincoming {};
	bool lforEgo {};
	bool lsqrtVersion {};
};

}

#endif /* COVARIATEDEGREEFUNCTION_H_ */
