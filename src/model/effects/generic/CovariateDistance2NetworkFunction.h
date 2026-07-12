/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2NetworkFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2NetworkFunction class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2NETWORKFUNCTION_H_
#define COVARIATEDISTANCE2NETWORKFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"

namespace siena
{

class CovariateDistance2NetworkFunction: public CovariateNetworkAlterFunction
{
public:
	CovariateDistance2NetworkFunction(std::string networkName, std::string covariateName,
						bool excludeMissing, bool outgoing, bool nc = false);
	virtual ~CovariateDistance2NetworkFunction();
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);


protected:
	bool missingDummy(int i) const;
	double averageAlterValue(int i) const;
	double totalAlterValue(int i) const;
	bool missingInDummy(int i) const;
	double averageInAlterValue(int i) const;
	double totalInAlterValue(int i) const;
	double similarityAvAlt(int i, int j) const;
	double varOutAvSimilarity(int i, int j) const;
	double varInAvSimilarity(int i, int j) const;

	// Whether _nc variants use raw (non-centered) covariate values.
	// Single source of truth for the distance-2 family; read by derived
	// alter/in-alter functions.
	inline bool nc() const;

private:
	bool lexcludeMissing {};
	bool loutgoing {};
	bool lnc {};
	double * laverageAlterValues {};
	double * ltotalAlterValues {};
	bool * laverageAlterMissing {};
	double * laverageInAlterValues {};
	double * ltotalInAlterValues {};
	bool * laverageInAlterMissing {};

};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns whether _nc variants use raw (non-centered) covariate values.
 */
bool CovariateDistance2NetworkFunction::nc() const
{
	return this->lnc;
}

}

#endif /* COVARIATEDISTANCE2NETWORKFUNCTION_H_ */
