/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2InAlterNetworkFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2InAlterNetworkFunction class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2INALTERNETWORKFUNCTION_H_
#define COVARIATEDISTANCE2INALTERNETWORKFUNCTION_H_

#include "CovariateDistance2NetworkFunction.h"

namespace siena
{

class CovariateDistance2InAlterNetworkFunction: public
	CovariateDistance2NetworkFunction
{
public:
	CovariateDistance2InAlterNetworkFunction(std::string networkName,
		std::string covariateName, bool excludeMissing, bool total);
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);
	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
	bool ltotal {};
};

}

#endif /* COVARIATEDISTANCE2INALTERNETWORKFUNCTION_H_ */
