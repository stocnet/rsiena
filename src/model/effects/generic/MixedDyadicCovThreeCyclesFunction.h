/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedDyadicCovThreeCyclesFunction.h
 *
 * Description: This file contains the definition of the
 * MixedDyadicCovThreeCyclesFunction class.
 *****************************************************************************/

#ifndef MIXEDDYADICCOVTHREECYCLESFUNCTION_H_
#define MIXEDDYADICCOVTHREECYCLESFUNCTION_H_

#include "DyadicCovariateMixedNetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: MixedDyadicCovThreeCyclesFunction class
// ----------------------------------------------------------------------------

/**
 * For sharedTo effect (see manual).
 */
class MixedDyadicCovThreeCyclesFunction : public DyadicCovariateMixedNetworkAlterFunction
{
public:
	MixedDyadicCovThreeCyclesFunction(std::string firstNetworkName,
			std::string secondNetworkName, 
			std::string dyadicCovariateName, double parameter, bool average);
			
	virtual ~MixedDyadicCovThreeCyclesFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
	virtual void preprocessEgo(int ego);
	virtual double value(int alter) const;

private:
	bool lroot {}; // should the square root be taken?
	bool lFirstWeight {}; // should there be centering?
	bool laverage {}; // should there be division by number of in-two-stars in first network?
	bool lSecondWeight {}; // average observed number of in-two-stars in first network
	std::string lvariableName {}; // name of first network
	SqrtTable * lsqrtTable;
	// for use in preprocessing:
	double * ltimesFound {};
	double lsumDegs {};
	int ln {};
};

}

#endif /* MIXEDDYADICCOVTHREECYCLESFUNCTION_H_ */
