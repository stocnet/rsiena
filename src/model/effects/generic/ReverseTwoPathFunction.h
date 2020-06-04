/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReverseTwoPathFunction.h
 *
 * Description: This file contains the definition of the
 * ReverseTwoPathFunction class.
 *****************************************************************************/

#ifndef REVERSETWOPATHFUNCTION_H_
#define REVERSETWOPATHFUNCTION_H_

#include "OneModeNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

class ConfigurationTable;


/**
 * Defines a function that returns the number of two-paths between
 * the ego and alters in a network of the given name.
 */
class ReverseTwoPathFunction:
	public OneModeNetworkAlterFunction, IntAlterFunction
{
public:
	ReverseTwoPathFunction(std::string networkName);

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

#endif /* REVERSETWOPATHFUNCTION_H_ */
