/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwCovariateAlterFunction.h
 *
 * Description: This file contains the definition of the
 * GwCovariateAlterFunction class.
 *****************************************************************************/

#ifndef GWBEHAVIORALTERFUNCTION_H_
#define GWBEHAVIORALTERFUNCTION_H_

#include "CovariateDistance2NetworkFunction.h"

namespace siena
{

class GwCovariateAlterFunction : public CovariateDistance2NetworkFunction
{
public:
	GwCovariateAlterFunction(std::string networkName,
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

#endif /* GWBEHAVIORALTERFUNCTION_H_ */
