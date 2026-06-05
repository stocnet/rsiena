/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwCovariateInAlterFunction.h
 *
 * Description: This file contains the definition of the
 * GwCovariateInAlterFunction class.
 *****************************************************************************/

#ifndef GWBEHAVIORINALTERFUNCTION_H_
#define GWBEHAVIORINALTERFUNCTION_H_

#include "CovariateDistance2NetworkFunction.h"

namespace siena
{

class GwCovariateInAlterFunction : public CovariateDistance2NetworkFunction
{
public:
	GwCovariateInAlterFunction(std::string networkName,
		std::string covariateName, bool excludeMissing, double parameter);

	virtual double value(int alter) const;

private:
	double gwWeight(double total) const;

	bool lexcludeMissing {};
	double lweight {};
	double lexpmweight {};
	double lexpfactor {};
};

}

#endif /* GWBEHAVIORINALTERFUNCTION_H_ */
