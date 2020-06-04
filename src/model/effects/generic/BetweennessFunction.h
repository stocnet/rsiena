/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BetweennessFunction.h
 *
 * Description: This file contains the definition of the
 * BetweennessFunction class.
 *****************************************************************************/

#ifndef BETWEENNESSFUNCTION_H_
#define BETWEENNESSFUNCTION_H_

#include "OneModeNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

class ConfigurationTable;


/**
 * Defines a function that returns the number of non-transitive two-paths
 * in a network of the given name.
 */
class BetweennessFunction: public OneModeNetworkAlterFunction, IntAlterFunction
{
public:
	BetweennessFunction(std::string networkName);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter);
	virtual int intValue(int alter);

private:
	ConfigurationTable * lpTable;
};

}

#endif /* BETWEENNESSFUNCTION_H_ */
