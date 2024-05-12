/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateMixedTwoPathFunction.h
 *
 * Description: This file contains the definition of the
 * SameCovariateMixedTwoPathFunction class.
 *****************************************************************************/

#ifndef SAMECOVARIATEMIXEDTWOPATHFUNCTION_H_
#define SAMECOVARIATEMIXEDTWOPATHFUNCTION_H_

#include "CovariateMixedNetworkAlterFunction.h"

namespace siena
{

class SameCovariateMixedTwoPathFunction: public CovariateMixedNetworkAlterFunction
{
public:
	SameCovariateMixedTwoPathFunction(std::string firstNetworkName,
		std::string secondNetworkName, std::string covariateName, 
		bool same, bool excludeMissing);
		
	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lsame {};
	bool lexcludeMissing {};
};

}

#endif /* SAMECOVARIATEMIXEDTWOPATHFUNCTION_H_ */
