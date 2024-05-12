/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2SimilarityNetworkFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2SimilarityNetworkFunction class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2SIMILARITYNETWORKFUNCTION_H_
#define COVARIATEDISTANCE2SIMILARITYNETWORKFUNCTION_H_

#include "CovariateDistance2NetworkFunction.h"

namespace siena
{

class CovariateDistance2SimilarityNetworkFunction: public
	CovariateDistance2NetworkFunction
{
public:
	CovariateDistance2SimilarityNetworkFunction(std::string networkName,
		std::string covariateName, bool excludeMissing);

	virtual double value(int alter) const;

private:
	bool lexcludeMissing {};
};

}

#endif /* COVARIATEDISTANCE2SIMILARITYNETWORKFUNCTION_H_ */
