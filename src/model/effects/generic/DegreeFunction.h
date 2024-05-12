/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DegreeFunction.h
 *
 * Description: This file contains the definition of the DegreeFunction class.
 *****************************************************************************/

#ifndef DEGREEFUNCTION_H_
#define DEGREEFUNCTION_H_

#include "NetworkAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the average out-degrees of alters
 * in a network with the given name.
 */
class DegreeFunction: public NetworkAlterFunction
{
public:
	DegreeFunction(std::string networkName, double par);
	virtual double value(int alter) const;
	
private:
	double lp {};
};

}

#endif /* DEGREEFUNCTION_H_ */
