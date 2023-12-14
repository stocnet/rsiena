/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2AlterNetworkFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2AlterNetworkFunction class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2ALTERNETWORKFUNCTION_H_
#define COVARIATEDISTANCE2ALTERNETWORKFUNCTION_H_

#include "CovariateDistance2NetworkFunction.h"

namespace siena
{


class CovariateDistance2AlterNetworkFunction: public
	CovariateDistance2NetworkFunction
{
public:
	CovariateDistance2AlterNetworkFunction(std::string networkName,
		std::string covariateName, double parameter, bool excludeMissing, bool total);
	virtual ~CovariateDistance2AlterNetworkFunction();

	virtual double value(int alter) const;

private:
	double lparameter {};
	bool lexcludeMissing {};
	bool ltotal {};
};

}

#endif /* COVARIATEDISTANCE2ALTERNETWORKFUNCTION_H_ */
