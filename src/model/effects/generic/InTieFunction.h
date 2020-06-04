/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InTieFunction.h
 *
 * Description: This file contains the definition of the
 * InTieFunction class.
 *****************************************************************************/

#ifndef INTIEFUNCTION_H_
#define INTIEFUNCTION_H_

#include "OneModeNetworkAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the values of ties to the ego
 * from alters in a network with the given name.
 */
class InTieFunction: public OneModeNetworkAlterFunction
{
public:
	InTieFunction(std::string networkName);

	virtual double value(int alter);
};

}

#endif /* INTIEFUNCTION_H_ */
