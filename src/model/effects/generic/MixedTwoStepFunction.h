/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedTwoStepFunction.h
 *
 * Description: This file contains the definition of the
 * MixedTwoStepFunction class.
 *****************************************************************************/

#ifndef MIXEDTWOSTEPFUNCTION_H_
#define MIXEDTWOSTEPFUNCTION_H_

#include "MixedNetworkAlterFunction.h"
#include "IntAlterFunction.h"
#include "network/NetworkUtils.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;
class MixedConfigurationTable;

// ----------------------------------------------------------------------------
// Section: MixedTwoStepFunction class
// ----------------------------------------------------------------------------


/**
 * Defines a function that returns the number of two-paths between
 * the ego and alters in a network of the given name.
 */
class MixedTwoStepFunction:
public MixedNetworkAlterFunction, IntAlterFunction
{
public:
	MixedTwoStepFunction(std::string firstNetworkName, std::string secondNetworkName,
						Direction firstDirection, Direction secondDirection, double par);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter) const;
	virtual int intValue(int alter);

private:
	MixedConfigurationTable * lpTable;
	Direction ldirection1 {};
	Direction ldirection2 {};
	bool ltrunc {}; // should the value be truncated?
	bool lroot {}; // should the square root be taken?
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /* MIXEDTWOSTEPFUNCTION_H_ */
