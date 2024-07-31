/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InStarFunction.h
 *
 * Description: This file contains the definition of the
 * InStarFunction class.
 *****************************************************************************/


#ifndef INSTARFUNCTION_H_
#define INSTARFUNCTION_H_

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
// Section: InStarFunction class
// ----------------------------------------------------------------------------


/**
 * Defines a function that returns the number of in-stars between
 * the ego and alters in a network of the given name.
 */
class InStarFunction: public NetworkAlterFunction, IntAlterFunction
{
public:
	InStarFunction(std::string networkName, bool root);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter) const;
	virtual int intValue(int alter);

private:
	ConfigurationTable * lpTable;
	bool lroot {}; // should the square root be taken?
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /* INSTARFUNCTION_H_ */
