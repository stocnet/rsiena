/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InStarOnlyFunction.h
 *
 * Description: This file contains the definition of the
 * InStarOnlyFunction class.
 *****************************************************************************/


#ifndef INSTARONLYFUNCTION_H_
#define INSTARONLYFUNCTION_H_

#include "NetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;
class ConfigurationTable;

// ----------------------------------------------------------------------------
// Section: InStarOnlyFunction class
// ----------------------------------------------------------------------------


/**
 * Defines a function that returns the number of in-stars between
 * the ego and alters in a network of the given name.
 */
class InStarOnlyFunction: public NetworkAlterFunction, IntAlterFunction
{
public:
	InStarOnlyFunction(std::string networkName, bool root);
	InStarOnlyFunction(std::string networkName, bool root, const bool simulatedState);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
	virtual void initialize(const Data * pData,
		State * pState, State * pSimulatedState, int period, Cache * pCache);

	virtual double value(int alter) const;
	virtual int intValue(int alter);

private:
	ConfigurationTable * lpTable;
	bool lroot {}; // should the square root be taken?
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /* INSTARONLYFUNCTION_H_ */
