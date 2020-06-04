/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutDegreeFunction.h
 *
 * Description: This file contains the definition of the
 * OutDegreeFunction class.
 *****************************************************************************/

#ifndef OUTDEGREEFUNCTION_H_
#define OUTDEGREEFUNCTION_H_

#include "NetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the out-degrees of alters in a network
 * with the given name.
 */
class OutDegreeFunction: public NetworkAlterFunction, IntAlterFunction
{
public:
	OutDegreeFunction(std::string networkName);

	virtual double value(int alter);
	virtual int intValue(int alter);
};

}

#endif /* OUTDEGREEFUNCTION_H_ */
