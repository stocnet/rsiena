/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoInDegreeFunction.h
 *
 * Description: This file contains the definition of the
 * EgoInDegreeFunction class.
 *****************************************************************************/

#ifndef EGOINDEGREEFUNCTION_H_
#define EGOINDEGREEFUNCTION_H_

#include "OneModeNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the in-degree of the ego regardless
 * of the alter.
 */
class EgoInDegreeFunction: public OneModeNetworkAlterFunction, IntAlterFunction
{
public:
	EgoInDegreeFunction(std::string networkName);

	virtual double value(int alter);
	virtual int intValue(int alter);
};

}

#endif /* EGOINDEGREEFUNCTION_H_ */
