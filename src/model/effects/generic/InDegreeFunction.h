/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InDegreeFunction.h
 *
 * Description: This file contains the definition of the
 * InDegreeFunction class.
 *****************************************************************************/

#ifndef INDEGREEFUNCTION_H_
#define INDEGREEFUNCTION_H_

#include "NetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the in-degrees of alters in a network with
 * the given name.
 */
class InDegreeFunction: public NetworkAlterFunction, IntAlterFunction
{
public:
	InDegreeFunction(std::string networkName);

	virtual double value(int alter);
	virtual int intValue(int alter);
};

}

#endif /* INDEGREEFUNCTION_H_ */
