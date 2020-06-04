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

class ConfigurationTable;


/**
 * Defines a function that returns the number of in-stars between
 * the ego and alters in a network of the given name.
 */
class InStarFunction: public NetworkAlterFunction, IntAlterFunction
{
public:
	InStarFunction(std::string networkName);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter);
	virtual int intValue(int alter);

private:
	ConfigurationTable * lpTable;
};

}

#endif /* INSTARFUNCTION_H_ */
