/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedTwoPathFunction.h
 *
 * Description: This file contains the definition of the
 * MixedTwoPathFunction class.
 *****************************************************************************/

#ifndef MIXEDTWOPATHFUNCTION_H_
#define MIXEDTWOPATHFUNCTION_H_

#include "MixedNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

class MixedConfigurationTable;


/**
 * Defines a function that returns the number of two-paths between
 * the ego and alters in a network of the given name.
 */
class MixedTwoPathFunction:
public MixedNetworkAlterFunction, IntAlterFunction
{
public:
	MixedTwoPathFunction(std::string firstNetworkName, std::string secondNetworkName);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter);
	virtual int intValue(int alter);

private:
	MixedConfigurationTable * lpTable;
};

}

#endif /* MIXEDTWOPATHFUNCTION_H_ */
