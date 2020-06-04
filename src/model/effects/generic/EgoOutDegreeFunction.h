/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoOutDegreeFunction.h
 *
 * Description: This file contains the definition of the
 * EgoOutDegreeFunction class.
 *****************************************************************************/

#ifndef EGOOUTDEGREEFUNCTION_H_
#define EGOOUTDEGREEFUNCTION_H_

#include "NetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the out-degree of the ego regardless
 * of the alter.
 */
class EgoOutDegreeFunction: public NetworkAlterFunction, IntAlterFunction
{
public:
	EgoOutDegreeFunction(std::string networkName);

	virtual double value(int alter);
	virtual int intValue(int alter);
};

}

#endif /* EGOOUTDEGREEFUNCTION_H_ */
