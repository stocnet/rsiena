/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HomCovariateMixedTwoPathFunction.h
 *
 * Description: This file contains the definition of the
 * HomCovariateMixedTwoPathFunction class.
 * Contributed by Robert Hellpap.
 *****************************************************************************/

#ifndef HOMCOVARIATEMIXEDTWOPATHFUNCTION_H_
#define HOMCOVARIATEMIXEDTWOPATHFUNCTION_H_

#include "CovariateMixedNetworkAlterFunction.h"

namespace siena
{

class HomCovariateMixedTwoPathFunction: public CovariateMixedNetworkAlterFunction
{
public:
	HomCovariateMixedTwoPathFunction(std::string firstNetworkName,
			std::string secondNetworkName, std::string covariateName,
			bool excludeMissing);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
};

}

#endif /* HOMCOVARIATEMIXEDTWOPATHFUNCTION_H_ */
