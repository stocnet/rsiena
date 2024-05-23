/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedThreeCyclesFunction.h
 *
 * Description: This file contains the definition of the
 * ThreeCyclesFunction class.
 *****************************************************************************/

#ifndef MIXEDTHREECYCLESFUNCTION_H_
#define MIXEDTHREECYCLESFUNCTION_H_

#include "MixedNetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;
class ConfigurationTable;

// ----------------------------------------------------------------------------
// Section: MixedThreeCyclesFunction class
// ----------------------------------------------------------------------------

/**
 * For sharedTo effect (see manual).
 */
class MixedThreeCyclesFunction : public MixedNetworkAlterFunction
{
public:
	MixedThreeCyclesFunction(std::string firstNetworkName,
			std::string secondNetworkName, double parameter, bool average);
	virtual ~MixedThreeCyclesFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter) const;

private:
	bool lroot {}; // should the square root be taken?
	bool lcenter {}; // should there be centering?
	bool laverage {}; // should there be division by number of in-two-stars in first network?
	double lavInTwoStar {}; // average observed number of in-two-stars in first network
	std::string lvariableName {}; // name of first network
	ConfigurationTable * lpFirstInStarTable;
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
	// for use in preprocessing:
	int * ltimesFound {};
	double lsumDegs {};
	int ln {};
};

}

#endif /*MIXEDTHREECYCLESFUNCTION_H_*/
